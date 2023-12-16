# SolarKAT
MeerKAT as a solar telescope.  

Radio frequency interference is a growing problem in radio astronomy, especially for new-generation telescopes such as MeerKAT. Its wide field of view and high sensitivity make the MeerKAT telescope capable of capturing the Sun even when looking far away from it (close to 90 degrees). This ability facilitates also contamination from the out-of-field Sun, causing, in some conditions, data corruption (which we call solar radio interference). To address this issue we present the SolarKAT pipeline. SolarKAT is a solar imaging pipeline that aims to mitigate solar interference in MeerKAT data.


## Installation from the repository
 You can get SolarKAT by cloning it from the repository to your local computer. To Run SolarKAT, you must first install Stimela (the framework on which SolarKAT is based). Stimela is installed in a virtual environment.
```
git clone https://github.com/ratt-ru/solarkat.git
```
## Running SolarKAT 

Before running SolarKAT, activate the Stimela virtual environment. Then run:
```
stimela run recipe.yml [recipe_name] obs=obs
```
Example: 
```
stimela run solarkat.yaml solarkat obs=L1
```

More details about the  SolarKAT pipeline can be found in the documentation here https://solarkat-docs.readthedocs.io/en/latest/index.html.
