# SolarKAT

SolarKAT is a solar imaging pipeline developed to mitigate solar interference in MeerKAT data and recover the visibilities rather than discarding them.

Radio frequency interference is a growing problem in radio astronomy, especially for new-generation telescopes such as MeerKAT. Its wide field of view and high sensitivity make the MeerKAT telescope capable of capturing the Sun even when looking far away from it. This ability facilitates contamination from the out-of-field Sun, causing, in some conditions, data corruption.

## Installation 

### Git clone
You can have access to SolarKAT by cloning it from the repository to your local computer. To Run SolarKAT, you must first install Stimela (the framework on which SolarKAT is based). Stimela is installed in a virtual environment.

```
git clone https://github.com/ratt-ru/solarkat.git

```

## Running SolarKAT 

Before running SolarKAT, create and activate the Stimela virtual environment. 

```
$: virtualenv -p python3 stimela_env
source stimela-env/bin/activate
```
Then run:
```
stimela run recipe.yml [recipe_name] obs=obs
```
Example: 
```
stimela run solarkat.yaml solarkat obs=L1
```

More details can be found in the Documentation here https://solarkat-docs.readthedocs.io/en/latest/index.html.


## Citing this Pipeline

We kindly request that any work using this pipeline cites the following reference:
```
@MISC{samboco2024solarkat,
    author = {Samboco, Victória G. and  Heywood, Ian and Smirnov, Oleg},
    title = {SolarKAT: Solar Imaging Pipeline for MeerKAT},
    keywords={Software},
    howpublished = {Astrophysics Source Code Library, record ascl:2401.013},
    eid = {ascl:2401.013}
    pages = {ascl:2401.013},
    archivePrefix = {ascl},
    eprint = {2401.013},
    month={January},
    year = {2024},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2024ascl.soft01013S}
    }
```
Thank you for considering our request.


## Note

This package is in public Beta stage. SolarKAT is feature-complete, has undergone testing on multiple datasets, and is suitable for broader use. While efforts have been made to ensure stability, there may still be some undiscovered issues. Users are encouraged to use it, provide feedback, and report any issues on the project's GitHub page (https://github.com/ratt-ru/solarkat).
