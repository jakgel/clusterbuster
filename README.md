# [ClusterBuster]()

A python 2.7 environment to jointly analyse the data products from cosmological simulations and radio sky surveys. It provides the functionalities to
go from data-cubes (simulations) or .fits maps (surveys) to object catalogues.

This python implementation provides a joint analysis framwork for radio relic surveys and simulations. You have to supply the map files or simulations cubes by yourself.
Because of the widely varying output of different simulations you also have write your own wrapper to parse simulaion output into this tool. ClusterBuster provides some functionalities for this purpose.
Once you have done this you are able to create a survey catalogue (.pickle) that you can use for further analysis. Approximate Bayesian computation for your own specified model is implemented in an submodel.
``astroABC`` requires ``NumPy``, and ``astropy``. ``astroABC``,``abcpmc``, and ``multiprocessing`` are optional and usefull if you want to analyse simulations.

In the subdirectories several functionalities are added.


##  Likely changes
- Remove cosmocalc; instead include it in the installing descriptions
- Change files paths in some .py with personal name folders, make them more generall (example RelicSurveys/RelicExtraction.py)
- Remove some DEBUGGING parts within the code
- Add the cosmic web datacubes examples and documentation including a link to the download directory
- Remove the NVSS images (currently 50MB) and add a description where to find them

Input for simulations:
- Cosmic web (cite, also cite the data science centre because it has the radio emission)
- MUSIC-2: threehundred massive cluster (needs processed shocks files which are not publicy available. Ask authors of [Nuza et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.470..240N))

Radio Surveys:
- NVSS (cite paper adsabs + add link to the personal catalogue of point sources+relic region+cluster information)


*To add:*
In the provided examples the focus on so called diffuse radio emision, which is synchrotron emission emerging from cosmic rays in galaxy clusters.


**[Disclaimer](#disclaimer)** |
**[Installation](#documentation)** |
**[Documentation](#documentation)** |
**[License](#license)** |

## Disclaimer
This .git outlines the overall strcuture of the software. As no larger release is planned the author currently regards this as an open documentation rather than an
fully supported openly accesible software. This means, that the software published here requires some cleanups of the code.


## Installation
Install via pip:

    pip install clusterbuster

Alternatively clone/fork this directory from  github, 
Alternatively also provide setup.sh?

## Dependencies
-NumPy >=1.8
-scipy >=0.16
-aplpy >=1.1
-NFW   >=0.2
-astropy >=1.3

*Additional packaes needed?*
Install addition packages via
conda install ephem
conda install -c astropy pyregion=1.2 
cosmocalc? --> change to astropy cosmology!



## Documentation
Currently sparse within the code


## Licence
Copyright 2018 Jakob Gelszinnis

ClusterBuster is free software made available under the MIT License. For details see the LICENSE.txt file.
