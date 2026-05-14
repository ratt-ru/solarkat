FROM ubuntu:24.04

LABEL org.opencontainers.image.source="https://github.com/ratt-ru/solarkat"
LABEL org.opencontainers.image.description="SolarKAT - Solar imaging pipeline for MeerKAT"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# System dependencies:
#   wsclean                         - imager used by stimela cabs (Ubuntu .deb omits chgcentre)
#   casacore-dev / casacore-data    - casacore headers + measures tables (for chgcentre)
#   libgsl-dev / liblapack-dev      - chgcentre link-time deps
#   build-essential                 - g++ to compile chgcentre
#   python3 / python3-venv          - Python 3.12 satisfies the >=3.11 requirement
#   git                             - chgcentre source pull + msutils pip git dep
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        wsclean \
        casacore-dev \
        casacore-data \
        libgsl-dev \
        liblapack-dev \
        build-essential \
        python3 \
        python3-venv \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Build chgcentre directly from wsclean v3.4 sources.
# Ubuntu's wsclean .deb (3.4-3) does not ship the chgcentre binary even though
# the upstream CMake installs it. We compile only the two chgcentre TUs against
# system casacore/GSL/LAPACK and the header-only aocommon submodule, avoiding
# the full wsclean CMake configure (HDF5/FFTW/schaapcommon/radler/...).
ARG WSCLEAN_REF=v3.4
RUN git clone --depth 1 --branch "${WSCLEAN_REF}" \
        https://gitlab.com/aroffringa/wsclean.git /tmp/wsclean \
    && git -C /tmp/wsclean submodule update --init --depth 1 external/aocommon \
    && g++ -std=c++17 -O2 \
        -I/tmp/wsclean/external/aocommon/include \
        /tmp/wsclean/chgcentre/main.cpp \
        /tmp/wsclean/chgcentre/progressbar.cpp \
        -lcasa_ms -lcasa_measures -lcasa_tables -lcasa_casa \
        -lgsl -lgslcblas -llapack \
        -o /usr/local/bin/chgcentre \
    && rm -rf /tmp/wsclean

# Sanity check
RUN command -v chgcentre

# uv for fast package installation
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Isolated venv (avoids PEP 668 friction with the system Python on Ubuntu 24.04)
RUN uv venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH" \
    VIRTUAL_ENV=/opt/venv

# Package source
COPY pyproject.toml README.md ./
COPY solarkat/ solarkat/

# Install solarkat plus its declared dependencies
# (stimela, cult-cargo, python-casacore, astropy, msutils)
RUN uv pip install --no-cache .

CMD ["stimela", "--help"]
