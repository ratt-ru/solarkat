# Image for SolarKAT's own python steps (scan bookkeeping, Sun ephemeris, region
# files, copying the Sun model back into the MS). The radio-astronomy tools that
# the recipe calls (wsclean, chgcentre, CASA, QuartiCal, breizorro, msutils) run
# in their cult-cargo images, so they are not installed here.
#
# Build:  docker build -t ghcr.io/ratt-ru/solarkat:<version> .
# Published automatically by .github/workflows/docker.yml on release tags.
FROM python:3.12-slim

LABEL org.opencontainers.image.source="https://github.com/ratt-ru/solarkat"
LABEL org.opencontainers.image.description="SolarKAT - Solar imaging pipeline for MeerKAT"
LABEL org.opencontainers.image.licenses="MIT"

ENV PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /opt/solarkat
COPY pyproject.toml README.md ./
COPY solarkat/ solarkat/

# python-casacore wheels bundle casacore, so no system packages are needed.
RUN pip install --upgrade pip \
    && pip install . \
    && python -c "import solarkat, casacore.tables, astropy; print('solarkat OK, astropy', astropy.__version__)"

# Run from a neutral directory so the installed package is what gets imported.
WORKDIR /
CMD ["python", "-c", "import solarkat; print(solarkat.__file__)"]
