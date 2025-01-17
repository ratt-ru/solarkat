<div align="center">
 
  <img src="solarkat-logo/solarkat_logo.png" alt="SolarKAT Logo" width="550">
 
</div>


<h1>Welcome to SolarKAT 🚀</h1>

SolarKAT is a solar imaging pipeline developed to mitigate solar interference in MeerKAT data and recover the visibilities rather than discarding them.

Radio frequency interference is a growing problem in radio astronomy, especially for new-generation telescopes such as MeerKAT.
Its wide field of view and high sensitivity make the MeerKAT telescope capable of capturing the Sun even when looking far away from it.
 This ability facilitates contamination from the out-of-field Sun, causing, in some conditions, data corruption.

## Installation 

### Git clone
You can have access to SolarKAT by cloning it from the repository to your local computer.
To Run SolarKAT, you must first install Stimela (the framework on which SolarKAT is based). Stimela is installed in a virtual environment. In addition to stimela install cult-cargo. Both should be in the master branch.

### Create and activate a virtual environment
```
python3.9 -m venv solarkat-env
```
```
source solarkat-env/bin/activate

```
then install stimela and cult-cargo inside the virtualenv

```
pip install stimela
```

```
pip install cult-cargo
```
and clone the solarkat repository
```
git clone https://github.com/ratt-ru/solarkat.git
```


### Installation using pip 

If you want to use SolarKAT as a Python package in your own recipe, it can be installed using pip (the updated version is ```solarkat==1.0.4```):

```
pip install solarkat

```

## Running SolarKAT 

You can execute stimela by running the command:
```
stimela run recipe.yml [recipe_name] obs=obs
```
Example:
```
stimela run solarkat.yaml solarkat obs=L1
```

More details can be found in the Documentation here https://solarkat-docs.readthedocs.io/en/latest/index.html.

### Running solarkat on Ilifu
You can run solarkat on Ilifu by setting ```backend``` to ```slurm```. You can also set the ```srun_opts``` such as ```time``` and ```mem``` (memory to use). For these values, there's no default value.
```
   backend:
      slurm:
        enable: true
        srun_opts:
          time: 1-12
          mem: 200GB
          cpus-per-task: 1
``` 



## Citing this Pipeline

We kindly request that any work using this pipeline cites the following reference:
```
@article{samboco2024solarkat,
    author = {Samboco, Victória G. and  Heywood, Ian and Smirnov, Oleg},
    title = {SolarKAT: A Solar Imaging Pipeline for MeerKAT},
    journal = {Astrophysics Source Code Library},
    pages={ascl:2401.013},
    month={January},
    year = {2024}
    }
```
Thank you for considering our request.


## Note

This package is in public Beta stage. SolarKAT is feature-complete, has undergone testing on multiple datasets, and is suitable for broader use.
While efforts have been made to ensure stability, there may still be some undiscovered issues.
Users are encouraged to use it, provide feedback, and report any issues on the project's GitHub page (https://github.com/ratt-ru/solarkat).

