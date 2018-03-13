#!/usr/bin/env python

# based on modify_fits by Alexander Drabent
from __future__ import division,print_function

import pyfits
import numpy as np
import sys
import math
import pyutil_my.SURVEYutil.py      as SUut
from   cosmocalc import cosmocalc 
from   scipy import ndimage
import warnings

warnings.simplefilter('error')


def XRay_Analaysis_fits(fits, z, cent=[-1,-1]):
  
  cosmoPara   =  cosmocalc(z, H0=70, WM=0.29, WV=0.71)
  plsc        =  cosmoPara['PS_kpc']                    # kpc/''

  image  = pyfits.open(fits)['PRIMARY'].data
  header = pyfits.getheader(fits)

  spixel       = abs(header['CDELT1']*3600)   # ''/pixel
  aperture     = 500/cosmoPara['PS_kpc']      # ''   / 500kpc
  aperture_pix = aperture / spixel            # pixel/ 500kpc
  
  s_pixel      = [sPixel*plsc, sPixel]
  
  XRay_Analaysis_image(image, z, s_pixel, cent)
  
  
def XRay_Analaysis_image(IMG, z, s_pixel, cent=[-1,-1], Xbeam=30, steps = [1,2,3]):
  
  
  if 0 in steps:
    print('Not implemented yet')
    ### Finding the X-ray structure of interest
    ### Getting masked image and center of that structure
    
    
  if cent[0] == -1:
      cent =  SUut.brightestPixel(ndi.gaussian_filter(image, (0.5*Xbeam/s_pixel[1],0.5*Xbeam/s_pixel[1])))  
      print(cent)
  
  
  if 1 in steps:
    ### Calculation of concentration parameter
    ### assuming aperture is 500 kpc
    aperture_pix  = 500/s_pixel[0]
    image         = np.copy(IMG, order='K') # To prevent alteration on the original IMG
    shape         = image.shape
    

    r = aperture_pix
    y,x = np.ogrid[-cent[1] : shape[0] - cent[1], -cent[0] : shape[1] - cent[0]]
    mask = x*x + y*y  <= r*r
    sum_500 = np.sum(image[mask])

    r = aperture_pix / 5.
    y,x = np.ogrid[-cent[1] : shape[0] - cent[1], -cent[0] : shape[1] - cent[0]]
    mask = x*x + y*y  <= r*r
    sum_100 = np.sum(image[mask])

    c = sum_100 / sum_500

    print '###################################'
    print 'concentration parameter c =', c
    print 'radius of aperture:', aperture_pix, 'pixels'
    print '###################################'



  if 2 in steps:
  ### Calculation of centroid shift
  ### assuming aperture is 500 kpc
  
    r_list = np.arange(aperture_pix, 0.05 * aperture_pix, - 0.05 * aperture_pix)
    separation = []
    for i,r in enumerate(r_list):
            y,x = np.ogrid[-cent[1] : shape[0] - cent[1], -cent[0] : shape[1] - cent[0]]
            mask = x*x + y*y  <= r*r
            image[~mask] = 0
            sum_r = np.sum(image)
            centroid = ndimage.measurements.center_of_mass(image)           
            #centroid = np.zeros(2)
            #for x in np.arange(shape[1]):
                    #for y in np.arange(shape[0]):
                            #centroid[1] += image[y][x]*x / sum_r
                            #centroid[0] += image[y][x]*y / sum_r
                            #pass
                    #pass
            separation.append(math.sqrt((centroid[1] - cent[0])**2 + (centroid[0] - cent[1])**2))
            print 'Centroid', i, ':', centroid, '| sum_r:', sum_r, '| separation: ', separation[i], 'pixels'
            pass

    print '--> centroid shift w =', np.std(separation) / r_list[0]
    print '###################################'



  if 3 in steps:
  ### Calculation of power ratio
  ### assuming aperture is 500 kpc

  #print 'P_0:', sum_500 * math.log(aperture_pix)
    m = 3
    a = 0
    b = 0
    r = aperture_pix
    y,x = np.ogrid[-cent[1] : shape[0] - cent[1], -cent[0] : shape[1] - cent[0]]
    mask = x*x + y*y  <= r*r
    np.savetxt('mask.txt', mask)
    image[~mask] = 0
    for x in np.arange(shape[1]):
            for y in np.arange(shape[0]):
                    try:
                            R    = math.sqrt((x - cent[0]) ** 2 + (y - cent[1])**2)
                            phi  = math.acos((x - cent[0]) / R)
                            pass
                    except RuntimeWarning:
                            continue
                            pass
                    a += image[y][x] * (R ** m) * math.cos(m * phi)
                    b += image[y][x] * (R ** m) * math.sin(m * phi)
                    pass
            pass
          
    P_0 = (sum_500 * math.log(r))**2
    P_3 = 1/(2 * m**2 * r**(2 * m)) * (a**2 + b**2)
    print 'Power ratio P_3 / P_0 =', P_3 / P_0
    print '###################################'