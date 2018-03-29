# [ClusterBuster]()

A python 2.7 environment to jointly analyse the radio maps from cosmological simulations and radio sky surveys. It provides the functionalities to
go from data-cubes (simulations) or .fits images (surveys) to object catalogues.

This python implementation provides a joint analysis framwork for radio relic surveys and simulations. You have to supply the map files or simulations cubes by yourself.
Because of the widely varying output of different simulations you also have write your own wrapper to parse simulaion output into this tool. ClusterBuster provides some functionalities for this purpose.
Once you have done this you are able to create a survey catalogue (.pickle) that you can use for further analysis. Approximate Bayesian computation for your own specified model is implemented in an submodel.
``astroABC`` requires ``NumPy``, and ``astropy``. The modules ``astroABC``,``abcpmc``, and ``multiprocessing`` are optional and usefull if you want to analyse simulations.

In the subdirectories several functionalities are added. In the provided examples the focus is on so called diffuse radio emision, which is synchrotron emission emerging from cosmic rays in galaxy clusters.
A future release will add additional modules to analyse at least MUSIC-2.

Radio Surveys:
- NVSS ([Condon et al. 1998](http://adsabs.harvard.edu/abs/1998AJ....115.1693C) + add link to the TLS catalogue of point sources+relic region+cluster information)



##  Likely changes
- The [cosmocalc](http://cxc.harvard.edu/contrib/cosmocalc/) will be removed. Instead installing descriptions will be layed out.
- Files paths in some .py will be made more general (example RelicSurveys/RelicExtraction.py)
- DEBUGGING   parts within the code will be removed
- DEVELOPMENT parts within the code will be removed
- Remove the NVSS images (currently 50MB) and add a description where to find them

Input for simulations likely added:
- Cosmic web (cite, also cite the data science centre because it has the radio emission)
- MUSIC-2: three-hundred massive clusters. Ask authors of [Nuza et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.470..240N) for the processed files that include the infered shock strength and other hydrodynamical properties within the galaxy clusters.


**[Disclaimer](#disclaimer)** |
**[Installation](#documentation)** |
**[Documentation](#documentation)** |
**[License](#license)** |

## Disclaimer
This .git outlines the overall structure of the software. As no larger release is planned the author currently regards this as an open documentation rather than an fully supported openly accesible software. This means, that the software published here requires some changes outlined above to be fully user-friendly.


## Installation
Install via pip:

    pip install clusterbuster

Alternatively clone/fork this directory from  github. 

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

*Optional packages*
apcpmc
astroabc



## Documentation
Currently sparse within the files. You can get an overview of the idea by reading section 2 & 3 of [Nuza et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.470..240N). The main repository includes. *surveyexample* shows how to get all relevant information from radio relics in the NVSS survey.  *surveysim* (once made accesible) includes the code to extract relevant information from the MUSIC-2 survey once you have additional files. *inference* lets you apply some more modern inference mechanism of physical models on the code.



## Licence
Copyright 2018 Jakob Gelszinnis

ClusterBuster is free software made available under the MIT License. For details see the LICENSE.txt file.
