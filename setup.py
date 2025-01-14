from setuptools import setup, find_packages

# Read the long description from README.md
with open('README.md', 'r', encoding='utf-8') as file:
    long_description = file.read()

setup(
    name='solarkat',
    version='1.0.4',
    author='Victória da Graça G. Samboco',
    author_email='vicky.samboco@gmail.com',
    description='Solar imaging pipeline for solar interference mitigation in MeerKAT data.',
    url='https://github.com/ratt-ru/solarkat.git',
    project_urls={
        'Documentation': 'https://solarkat-docs.readthedocs.io/en/latest/',
        'Source Code': 'https://github.com/ratt-ru/solarkat',
    },
    packages=find_packages(),  
    long_description=long_description,
    long_description_content_type='text/markdown',

    package_data={
        'solarkat.solarkat-pipeline': [
            'solarkat.yaml',
            'solarkat-cabs.yaml',
            'solarkat-observation-sets.yaml',
        ],
    },

    py_modules=['solarkat.solarkat'],  

    install_requires=[
        'stimela>=2.0.1',
        'cult-cargo>=0.1.3',
        'quartical>=0.2.1',
        'msutils>=1.2.0',
    ],

    python_requires='>=3.9',  

    classifiers=[ 
        'Development Status :: 4 - Beta', 
        'Intended Audience :: Science/Research',  
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering :: Radioastronomy',
    ],

    include = [
    'solarkat/solarkat-pipeline/solarkat.yaml',
    'solarkat/solarkat-pipeline/solarkat-cabs.yaml',
    'solarkat/solarkat-pipeline/solarkat-observation-sets.yaml',
    'solarkat/solarkat/solarkat.py',
    ]
)

