#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:34:33 2017

I am quite certain that many functionalities are replaceble by modules that already exist.
The functionalitie sprovided here are mostly to modify numpy array maps, find contours, 

@author: jakobg

"""


from __future__ import division,print_function

import cv2
import warnings
import numpy             as np
import matplotlib.pyplot as plt 

from astropy.modeling import models, fitting


# Rotation angle in radians. The rotation angle increases !clockwise!, from the positive x-axis. Astropy 0.3
def ImageGaussian(shape, amp, stddev, center, theta=0, FWHM=False, order=2):

  if FWHM:
    stddev[:] = [x/1.17741/2 for x in stddev]  
    
  g = models.Gaussian2D(amplitude=amp, x_mean=center[1], y_mean=center[0], x_stddev=stddev[0], y_stddev=stddev[1], theta=np.radians(theta) )
  gauss  = np.zeros( (shape[0], shape[1])  ) 
    
  error  = 0
  range_xy    =  np.asarray(np.around(center , decimals = 0), int) 
  range_delta =  int(np.around(3.5*stddev[0], decimals = 0))
  for i in range(range_xy[0]-range_delta,range_xy[0]+range_delta):
     for j in range(range_xy[1]-range_delta,range_xy[1]+range_delta):
        try:
          gauss[j,i]  =  g(j,i)            #np.exp(  -np.power( (np.sqrt( (i-center[0])**2 + (j-center[1])**2  )) / (np.sqrt(2)*width)  , order )   )  
        except:
           error += 1
           
  gauss[np.isnan(gauss)] = 0         
  return  gauss

  
# Rotation angle in radians. The rotation angle increases counterclockwise, from the positive x-axis. Astropy 1.0
def ImageGaussian_inv(shape, amp, stddev, center, theta=0, FWHM=False, order=2):

  return  ImageGaussian(shape, amp, stddev, center, (450-theta)%180 , FWHM=FWHM, order=order) 
  

#Also scipy.optimize can provide good results, see http://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
def FitGaussian(img, amp, stddev, center, theta=0, FWHM=False): 

  if FWHM:
    stddev[:] = [x/1.17741/2 for x in stddev]  
    
  xy = np.where(img>0)
    
  # Fit the data using astropy.modeling
  p_init = models.Gaussian2D(amplitude=amp, x_mean=center[1], y_mean=center[0], x_stddev=stddev[0],y_stddev=stddev[1], theta=np.radians(theta) )
  fit_p  = fitting.LevMarLSQFitter()

  with warnings.catch_warnings():
      # Ignore model linearity warning from the fitter
      warnings.simplefilter('ignore')
      p = fit_p(p_init, xy[0], xy[1], img[xy])

  return p

   
def FitGaussian_inv(img, amp, stddev, center, theta=0, FWHM=False): 
  g = FitGaussian(img, amp, stddev, center, theta=theta, FWHM=FWHM)
  g.theta.value =  (2*np.pi-g.theta.value)%np.pi
  x = g.x_mean
  y = g.y_mean
  g.x_mean  = y
  g.y_mean = x
  return g

#====
from scipy.fftpack import fft, ifft
def FourierFilter(image):

  yf       = np.fft.ifft2(image)
  ffilter  = ImageGaussian( image.shape, 20, 2, [image.shape[0]/2-1,image.shape[1]/2-1] ) #  # 1 - ImageGaussian( image.shape, 20, 2, [image.shape[0]/2-1,image.shape[1]/2-1] )   
  filtered = np.multiply(yf,ffilter)
  
  
  return np.real(np.fft.fft2(filtered))


 
def is_contour_bad(c):
     # approximate the contour
     peri = np.arcLength(c, True)
     approx = np.approxPolyDP(c, 0.02 * peri, True)
 
     ## the contour is 'bad' if it is roughly an rectangle
     #return not len(approx) < 10 #  3.5
     
     # the contour is 'bad' if it is only one point
     return c.shape[0] == 0
     
     
     

def removingcontours(image):
 
  edged = np.Canny(image, 50, 100)
  np.imshow("Original", image)
 
  # find contours in the image and initialize the mask that will be
  # used to remove the bad contours
  (_, cnts, _) = np.findContours(edged.copy(), np.RETR_LIST, np.CHAIN_APPROX_SIMPLE)
  mask = np.ones(image.shape[:2], dtype="uint8") * 255
 
  # loop over the contours
  for c in cnts:
     # if the contour is bad, draw it on the mask
     if is_contour_bad(c):
        np.drawContours(mask, [c], -1, 0, -1)
 
  # remove the contours from the image and show the resulting images
  image = np.bitwise_and(image, image, mask=mask)
  #np.imshow("Mask", mask)
  #np.imshow("After", image)
  #np.waitKey(0)

  return image


def inverte(imagem, name):
    imagem = (255-imagem)
    np.imwrite(name, imagem)
    
    
def plot_comparison(original, filtered, filter_name):
  
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12, 4))
    ax1.imshow(original, cmap=plt.cm.gist_heat)
    ax1.set_title('original')
    ax1.axis('off')
    ax2.imshow(filtered, cmap=plt.cm.gist_heat)
    ax2.set_title(filter_name)
    ax2.axis('off')
    ax3.imshow(original-filtered, cmap=plt.cm.gist_heat)
    ax3.set_title('(relative) residuum')
    #ax3.axis('off')
    return fig

    
    
def ContourMasking(image, cnts):
  mask = np.zeros(image.shape[:2], dtype="uint8") 
 
  # loop over the contours
  for c in cnts:
        cv2.drawContours(mask, [c], -1, 1, -1)
  
  region = np.multiply(image,mask) # masking
  region[np.isnan(region)] = 0     # For contour masked  NVSS images I encountered the isue that some values where nan

  return region 
    

  
# mask image, so that every value outside the contour is set to zero
def GetThreshContours(image, thresh, contourMask=[]):
  import cv2 as cv2 #opencv as cv2
  
  if len(contourMask) > 0:
      image = ContourMasking(image, contourMask)
  
  Ndil   = 1

  #img2 =  np.greater(img, (img*0.)+cutoff) #img[img>cutoff]  mhhh, causes prblems
  #img3 =  morp.binary_dilation(img2, iterations=Ndil)
  
  #Deprecated: image = stats.threshold(image, threshmin=thresh, threshmax=None)
  image[image < thresh]    = 0
#  thresholded              = image[image > thresh] # image.astype(np.uint8) 
#  image[image > thresh]    = 1         

        
  'Tresholded is the image in mutliples of treshold. This step is quite costly, keeping in mind that we just want to get the contours'      
  thresholded             =  ( np.clip( image/thresh, 0, 1) ).astype(np.uint8)                                                          
  hierarchy, contours, _   = cv2.findContours(thresholded ,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
  
  return (contours, hierarchy)   



  
def PCA_analysis (image, upscale, cutoff, minArea, contourMask = [], Imcenter=[-1,-1]):
    ''' This was an experiment for the analysisi of X-ray maps. '''



    #=== Define Image parameters
    if Imcenter[0] == -1:
        Imcenter = [image.shape[0]/2,image.shape[1]/2]

        
    (contours, hierarchy)  = GetThreshContours(image, cutoff, upscale, contourMask)   
    while len(contours) == 0 or np.contourArea(sorted(contours, key=lambda cnt: np.contourArea(cnt), reverse=True)[0]) < minArea:
        print('No larger area above threshhold identified. Cutoff reduced by a factor of 0.8!')
        cutoff                = 0.8*cutoff
        (contours, hierarchy)  = GetThreshContours(image, cutoff, upscale, contourMask)   
    cnt = sorted(contours, key=lambda cnt: np.contourArea(cnt), reverse=True)[0]
    M = np.moments(cnt)

    if cnt.shape[0] < 2 :
        print('Contour of less than three point identified:', cnt, '. It will not be further considered.')
        return 
    if M['m00']    != 0 :
        cx = int(M['m10']/M['m00'])
        cy = int(M['m01']/M['m00'])

        # Miimal enclosing circle
        (x,y),radius = np.minEnclosingCircle(cnt)
        center = (x,y)
	
        # Fitting an ellipse
    if cnt.shape[0] > 4:
        ellipse = np.fitEllipse(cnt)
        a =  ellipse[1][1]
        b =  ellipse[1][0]
    else:
        a = 0
        b = 0

        # mask image, so that every value is set to zero
        mask    = np.zeros(image.shape[:2], dtype="uint8") 
        np.drawContours(mask, [cnt], -1, 1, -1)
        region  = np.multiply(image,mask)  
        region[np.isnan(region)] = 0     #For contour masked  NVSS images I encountered the isue that some values where nan
	  
        PCA_histo = np.where(region > cutoff/upscale )   #
        PCA_array = np.column_stack((PCA_histo[0],PCA_histo[1],image[PCA_histo]))    #np.hstack( (PCA_histo, image[PCA_histo]) )
        PCA_array_b = np.column_stack((PCA_histo[0],PCA_histo[1]))
          
  
        # compute flux and luminosity weighted coordinaetes/distance to cluster center
        flux       = np.sum(region)
      
        xmean    = np.average( PCA_array_b[:,0], weights = image[PCA_histo])
        ymean    = np.average( PCA_array_b[:,1], weights = image[PCA_histo])
        #mean   = np.mean(PCA_array, axis=0  )  

        plt.imshow(np.log10(region+0.2*image+1e-9))
        plt.show()


''' Based on the concept of Image moment (inckuding raw moemnts,  and central moments)
#===== from http://stackoverflow.com/questions/9005659/compute-eigenvectors-of-image-in-python/9007249#9007249
#=====    & http://stackoverflow.com/questions/5869891/how-to-calculate-the-axis-of-orientation
'''
def raw_moment(data, iord, jord):
    nrows, ncols = data.shape
    y, x = np.mgrid[:nrows, :ncols]
    data = data * x**iord * y**jord
    return data.sum()

    
def inertial_axis(data):
    """Calculate the x-mean, y-mean, and cov matrix of an image."""
    data_sum = data.sum()        
    m10 = raw_moment(data, 1, 0)
    m01 = raw_moment(data, 0, 1)
    x_bar = m10 / data_sum  # a normed raw moment
    y_bar = m01 / data_sum  # a normed raw moments
    u11 = (raw_moment(data, 1, 1) - x_bar * m01) / data_sum
    u20 = (raw_moment(data, 2, 0) - x_bar * m10) / data_sum
    u02 = (raw_moment(data, 0, 2) - y_bar * m01) / data_sum
    cov = np.array([[u20, u11], [u11, u02]])

    if (u20-u02) != 0:
        angle = 0.5 * np.arctan(2 * u11 / (u20 - u02))
    else:
        angle = 0
    if  u02 > u20:
      angle += np.pi/2
      
    return x_bar, y_bar, cov, angle
    

      
#def plot_subplot(paw, ax, center, cov):
#    ax.imshow(paw)
#    plot_bars(center[0], center[1], cov, ax) #plot_bars_simple(x_bar, y_bar, angle, ax) 
#
#    return

#from astropy.nddata import Cutout2D
#from astropy.nddata.utils import extract_array
#def plot(images, size=False, center=False):
#      
#    fig1     = plt.figure()
#    ax1      = fig1.add_subplot(1,1,1)
#    paw = images #.sum(axis=2)
#    x_bar, y_bar, cov, angle = inertial_axis(paw)
#    if not center:
#      center =   (y_bar, x_bar)
#    if size:       
#        #cutout     = Cutout2D(paw, center, size)
#        #paw   = cutout.data
#        #center     = cutout.position_cutout
#
#        paw  =  extract_array(paw, size, center)
#        center    =  ( int(paw.shape[1]/2),int(paw.shape[0]/2) )
#    plot_subplot(paw, ax1, center, cov)
#    
#
#    fig1.suptitle('Original')
    #fig2.suptitle('Rotated')

#def plot_bars(x_bar, y_bar, cov, ax):
#    """Plot bars with a length of 2 stddev along the principal axes."""
#    def make_lines(eigvals, eigvecs, mean, i):
#        """Make lines a length of 2 stddev."""
#        std = np.sqrt(eigvals[i])
#        vec = 2 * std * eigvecs[:,i] / np.hypot(*eigvecs[:,i])
#        x, y = np.vstack((mean-vec, mean, mean+vec)).T
#        return x, y
#    mean = np.array([x_bar, y_bar])
#    eigvals, eigvecs = np.linalg.eigh(cov)
#    ax.plot(*make_lines(eigvals, eigvecs, mean, 0), marker='o', color='white')
#    ax.plot(*make_lines(eigvals, eigvecs, mean, -1), marker='o', color='red')
#    ax.axis('image')
#    
#    
#def plot_bars_simple(x_bar, y_bar, angle, ax):
#    def plot_bar(r, x_bar, y_bar, angle, ax, pattern):
#        dx = r * np.cos(angle)
#        dy = r * np.sin(angle)
#        ax.plot([x_bar - dx, x_bar, x_bar + dx], 
#                [y_bar - dy, y_bar, y_bar + dy], pattern)
#    plot_bar(3, x_bar, y_bar, angle + np.radians(90), ax, 'wo-')
#    plot_bar(9, x_bar, y_bar, angle, ax, 'ro-')
#    ax.axis('image')


def normalize_image(image, smin=0, smax=1, valmax=0):
    a = np.min(image)
    #print np.max(image)
    
    if valmax != 0:
      b = valmax
    else:
      b = np.max(image)
      
    scaled = (image - a ) / ((b - a)) * (smax-smin) + smin
    scaled[np.where(scaled > smax)] = smax
    return scaled


def padwithtens(vector, pad_width, iaxis, kwargs):
    vector[:pad_width[0]]  = 0
    vector[-pad_width[1]:] = 0
    return vector