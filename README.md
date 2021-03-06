# [ClusterBuster]()

A python 2.7 environment to jointly analyse the radio maps from cosmological simulations and radio sky surveys. It provides the functionalities to
go from data-cubes (simulations) or .fits images (surveys) to object catalogues.

It provides a joint analysis framework for radio relic surveys and simulations. You have to supply the map files or simulations cubes by yourself.
Because of the widely varying output of different simulations you also have write your own wrapper to parse simulaion output into this tool. ClusterBuster provides some functionalities for this purpose.
Once you have done this you are able to create a survey catalogue (.pickle) that you can use for further analysis. Approximate Bayesian computation for your own specified model has to be implemented in a submodel (so some coding is required).
``astroABC`` requires ``NumPy``, and ``astropy``. The module [``abcpmc``](https://github.com/jakeret/abcpmc) is optional and usefull if you want to analyse simulations.

In the subdirectories several functionalities are added. In the provided examples the focus is on so-called diffuse radio emision, which is synchrotron emission emerging from cosmic rays in galaxy clusters.

##### Radio Surveys:
- NVSS ([Condon et al. 1998](http://adsabs.harvard.edu/abs/1998AJ....115.1693C) + add link to the TLS catalogue of point sources+relic region+cluster information)


##### Scripts for analysis of simulations is added:

- MUSIC-2: 283 massive clusters resimulated from Multidark. Ask authors of [Nuza et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.470..240N) for the processed files that include the infered shock strength and other hydrodynamical properties within the galaxy clusters.

####  Likely changes
- Installing descriptions will be layed out.
- NVSS poststamp images (currently 50MB) will be removed and a description will be added about  where to find them

**[Disclaimer](#disclaimer)** |
**[Installation](#documentation)** |
**[Documentation](#documentation)** |
**[License](#license)** |

## Disclaimer
This .git outlines the overall structure of the software. As no larger release is planned the author currently regards this as an open documentation rather than an fully supported openly accesible software. This means, that the software published here requires some changes outlined above to be fully user-friendly.


## Installation
Clone/fork this directory from  github. 
After download go to the directory and add the repository to your bashrc
    >> echo "export PYTHONPATH=\PYTHONPATH:$(pwd)" >> ~/.bashrc

## Dependencies
- NumPy >=1.8
- scipy >=0.16
- pandas >= 0.22
- opencv >= 3.4
- aplpy >=1.1
- NFW   >=0.2
- astropy >=1.3
- ephem >=3.7.6
- PyPDF2 >= 1.26.0
- matplotlib > (needs matplotlib.mlab)
- reproject (the NVSS subtraction part)
Older versions may work but are untested.


## Documentation
Currently sparse within the files. You can get an overview of the idea by reading section 2 & 3 of [Nuza et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.470..240N). The main repository includes. *surveyreal* shows how to get all relevant information from radio relics in the NVSS survey.  *surveysim* (once made accesible) includes the code to extract relevant information from the MUSIC-2 survey once you have additional files. *inference* lets you apply some more modern inference mechanism of physical models on the code.

![Example of parameter inference with [ABCPMC](https://github.com/jakeret/abcpmc)](inference/Example_abcpmc.png)

## Licence
ClusterBuster is free software made available under the MIT License. For details see the LICENSE.txt file.
