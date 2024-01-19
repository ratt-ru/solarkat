from setuptools import setup, find_packages

setup(
    name='solarkat',  # Replace with your own package's name
    version='0.1.0',  # Package version
    author='Victória Samboco',  # Your name
    author_email='vicky.samboco@gmail.com',  # Your email
    description='A solar imaging pipeline for solar interference mitigation in MeerkAT,  # Short description
    long_description=open('README.md').read(),  # Long description read from the README.md file
    long_description_content_type='text/markdown',  # Content type of the long description, typically text/markdown or text/x-rst
    url='https://github.com/ratt-ru/solarkat.git',  # URL to your package's repository
    packages=find_packages(),  # Automatically find all packages in the directory
    install_requires=[  # List of dependencies
        'aiobotocore==2.4.0',
        'aiohttp==3.8.1',
        'aioitertools==0.10.0',
        'aiosignal==1.2.0',
        'antlr4-python3-runtime==4.9.3',
        'appdirs==1.4.4',
        'asciitree==0.3.3',
        'astLib==0.11.7',
        'astro-kittens==1.4.6',
        'astro-tigger-lsm==1.7.2',
        'astropy==5.1',
        'async-timeout==4.0.2',
        'attrs==22.1.0',
        'bokeh==2.4.3',
        'botocore==1.27.59',
        'breizorro==0.1.0',
        'certifi==2022.9.24',
        'charset-normalizer==2.1.1',
        'cli-ui==0.17.2',
        'click==8.1.3',
        'cloudpickle==2.2.0',
        'codex-africanus==0.3.3',
        'colorama==0.4.5',
        'Columnar==1.4.1',
        'commonmark==0.9.1',
        'configparser==5.3.0',
        'contextlib2==21.6.0',
        'contourpy==1.0.5',
        '-e git+https://github.com/paoloserra/crystalball.git@1c24d1fbaea8a22c31a42a941103e305d9906d0f#egg=crystalball',
        'cycler==0.11.0',
        'dask==2022.9.2',
        'dask-ms==0.2.14',
        'decorator==5.1.1',
        'dill==0.3.5.1',
        'distributed==2022.9.2',
        'docopt==0.6.2',
        'donfig==0.7.0',
        'entrypoints==0.4',
        'fasteners==0.18',
        'fonttools==4.38.0',
        'frozenlist==1.3.1',
        'fsspec==2022.8.2',
        'future==0.18.2',
        'future-fstrings==1.2.0',
        'HeapDict==1.0.1',
        'idna==3.4',
        'importlib-metadata==5.0.0',
        'iniconfig==1.1.1',
        'Jinja2==3.1.2',
        'jmespath==1.0.1',
        'kiwisolver==1.4.4',
        'llvmlite==0.39.1',
        'locket==1.0.0',
        'loguru==0.6.0',
        'MarkupSafe==2.1.1',
        'matplotlib==3.6.1',
        'msgpack==1.0.4',
        'msutils==1.2.0',
        'multidict==6.0.2',
        'multiprocess==0.70.13',
        'munch==2.5.0',
        'nose==1.3.7',
        'numba==0.56.2',
        'numcodecs==0.10.2',
        'numpy==1.23.3',
        'omegaconf==2.2.3',
        'packaging==21.3',
        'pandas==1.5.0',
        'partd==1.3.0',
        'pathos==0.2.9',
        'Pillow==9.2.0',
        'pkg_resources==0.0.0',
        'pluggy==1.0.0',
        'pox==0.3.1',
        'ppft==1.7.6.5',
        'psutil==5.9.2',
        'py==1.11.0',
        'pydantic==1.10.2',
        'pyerfa==2.0.0.1',
        'Pygments==2.13.0',
        'pyparsing==3.0.9',
        'pytest==7.1.3',
        'python-casacore==3.5.2',
        'python-dateutil==2.8.2',
        'pytz==2022.4',
        'PyYAML==6.0',
        '-e git+https://github.com/ratt-ru/QuartiCal.git@353279b119cfbafcb53f16154bef41e16f254ca9#egg=quartical',
        'regions==0.5',
        'requests==2.28.1',
        'rich==12.5.1',
        'ruamel.yaml==0.17.21',
        'ruamel.yaml.clib==0.2.6',
        's3fs==2022.8.2',
        'scabha==2.0.0',
        'schema==0.7.5',
        'scipy==1.9.1',
        'six==1.16.0',
        'sortedcontainers==2.4.0',
        '-e git+https://github.com/caracal-pipeline/stimela2.git@3a631db55ee3c9d758c7f55300fb8de7716d38de#egg=stimela',
        'tabulate==0.8.10',
        'tblib==1.7.0',
        'tbump==6.9.0',
        'tomli==2.0.1',
        'tomlkit==0.11.5',
        'toolz==0.12.0',
        'tornado==6.1',
        'typing_extensions==4.3.0',
        'Unidecode==1.3.6',
        'urllib3==1.26.12',
        'wcwidth==0.2.5',
        'wrapt==1.14.1',
        'xarray==2022.9.0',
        'yarl==1.8.1',
        'zarr==2.13.2',
        'zict==2.2.0',
        'zipp==3.8.1',
        '-e git+https://github.com/ratt-ru/breizorro.git#egg=breizoirro'

      # Replace these with your actual dependencies
      # Add more dependencies as needed
    ],
    classifiers=[  # Classifiers help categorize your project and give it more visibility
        'Development Status :: 3 - Alpha',  # Choose the appropriate development status
        'Intended Audience :: Radio astronomers, Developers',  # Define the audience for your package
        'License :: OSI Approved :: MIT License',  # Choose the appropriate license
        'Programming Language :: Python :: 3',  # Supported Python versions
        'Programming Language :: Python :: 3.8',  # Specify specific Python versions if needed
        # Add more classifiers as appropriate for your package
    ],
    python_requires='>=3.8',  # Minimum version requirement of the Python for your package
    entry_points={  # Optional: Define console scripts or GUI scripts here
        'console_scripts': [
            'your_command = your_package.module:function',
            # Example: 'pip = pip._internal.cli.main:main'
        ]
    },
)
