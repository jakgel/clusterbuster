#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 17:03:47 2017

@author: jakobg
"""


#  Units, conversion factors, and input variables
radiounit_A  =  1                     # W/Hz        --- Unit of particle luminousity in .radio snaps, formerly 1e23 for Sebastians cubes
Jy           =  1e-23                 # erg/s/Hz/cm^2
Jy_SI        =  1e-26                 # W/Hz/m^2
H0           =  70                    # 0.7  ... Hubble constant used in MUSIC-2
FWHM2sigma   =  1/2.354               # convert FWHM to sigma for a gaussian