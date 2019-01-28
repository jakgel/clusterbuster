#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 10:52:22 2017

@author: jakobg
"""
from __future__ import division,print_function

import glob
import os                                   
import numpy                   as np
import clusterbuster.iout.misc as iom

def AddFilesToSurvey(survey, savefolder, verbose = True, clusterwise=False):
    """
    Adds galaxy cluster in a specified folder to an survey.
    
    Currently clusterwise=True is the case for the normal Run() and clusterwise=False is the case for the ABC-routine.
    How about CW?
    """

    """ This replaces tha galaxy clusters of a survey object with all pickled galaxy clusters in an particular 'savefolder' """
    minsize  = 1
    location = savefolder+'/pickled/'
    location = location.replace('//', '/')
    GCls     = []

    if clusterwise: 
        fn = 'GCl'
    else:
        fn = 'MonsterPickle'
        

    if verbose: print('%s%s*.pickle' % (location, fn))
        
    for filename in glob.glob('%s%s*.pickle' % (location,fn)):  #glob.glob('%srelics/*.pickle' % (location)):
        if verbose:
            print('surveyutil::AddFilesToSurvey()::', filename)
        if os.path.getsize(filename) > minsize:
             if verbose: print('surveyutil::AddFilesToSurvey()::filename', filename)
             items = iom.unpickleObjectS(filename)
             for (Cluster, Rmodel) in items:
                 GCls.append(Cluster)
             os.remove(filename) 

    GCls = sorted(GCls, key= iom.Object_natural_keys)
    survey.GCls = GCls
     
    if len(GCls) == 0:
        print('surveyutil::AddFilesToSurvey() len(GCls)', len(GCls))
    
    iom.pickleObject(survey, location, 'Survey')
    return survey

def TestPar(var):
  
  return var in ['True', '1', 'TRUE', 'true', True]



def interpret_parset(parfile, repository='/parsets/', default='default.parset', verbose=False, relative=None, oldstyle=False):
#This set loads a parset and replaxes one with the default one, alternatively you can give the default one and then modify the entries by another parset2dict
    from pathlib2 import Path
    import os


    if relative is None:
        if  Path(os.getcwd() + repository + default).is_file():
            relative = True
        elif Path(repository + default).is_file():
            relative = False
             
    def_dict  = iom.parset2dict(repository + default, relative=relative) # relative=relative
    new_dict  = iom.parset2dict(repository + parfile, relative=relative)
     
    
    comb_dict = def_dict.copy()
    comb_dict.update(new_dict)
    
    if verbose: 
        print(parfile, comb_dict)


    z_arr    =  [float(z) for z in   iom.str2list(comb_dict['z_range']) ]
    if oldstyle:
        # Now this parset is used to create some lists needed for this script
        B0_arr   =  iom.createlist(iom.str2list(comb_dict['B0']       ),comb_dict['B0_N']     , interpol= comb_dict['B0_type'])
        nu_arr   =  iom.createlist(iom.str2list(comb_dict['nu']       ),comb_dict['nu_N']     , interpol= comb_dict['nu_type'])
        eff_arr  =  iom.createlist(iom.str2list(comb_dict['eff_range']),comb_dict['eff_steps'], interpol= comb_dict['eff_type'])[::-1]  #invert array
        returnargs = (comb_dict, B0_arr, nu_arr, eff_arr, z_arr)
    else:
        returnargs = (comb_dict, z_arr)



    return returnargs
     
""" former CW """


def brightestPixel(image):  # Runs into issues if there is more than one brightest pixel, this is also why a smoothed image is prefered

    VMax = image.max() 
    return np.where(image == VMax)



def cmask(image,radius,center=None):
    """ takes the 2D-array and creates a mask with a radius of 'radius' elements
    If not specified by index the center of the image is given as the center.
    inspired by http://stackoverflow.com/questions/8647024/how-to-apply-a-disc-shaped-mask-to-a-numpy-array 
    Returns: ... """
    

    
    if center is None:
        center = int(image.shape[0]/2.), int(image.shape[1]/2.)


    a, b = center[0], center[1]
    n1 = image.shape[0]
    n2 = image.shape[1]

    
    y,x = np.ogrid[-a:n1-a, -b:n2-b]
    mask = x*x + y*y >= radius*radius

    return mask  


def gaussian_pseudo(x,x0,sigma):
# the integral is not zero, however we only compare absolute values which is fine
  return  np.exp(-np.power((x - x0)/sigma, 2.)/2.)  #1./(np.sqrt(2.*np.pi)*sigma)*

def Recenter(GCl, image, subp, image2mask=None, setto=0.): # a function of galaxy cluster?
    """ recenters the cluster centric coordinates basede on the mass proxy 
        WORK IN PROGRESS: Because It changes the GCl object it should be one of it's functionalities!
        """
    # draw region of Mwir around current center """
    # find maximum within Rvir ...
    radius = GCl.R200/(146.48/subp)
    mask      = cmask(image,radius)  #(1+GCl.z)
    masked_im = np.copy(image)
    masked_im[mask] = -10
    cen_x, cen_y = brightestPixel(masked_im)
    off_x = cen_x[0] - int(masked_im.shape[0]/2)
    off_y = cen_y[0] - int(masked_im.shape[1]/2)
    
    GCl.RA.value       = GCl.RA  - off_y* 146.48/(subp*GCl.cosmoPS*3600)/(1+GCl.mockobs.z_snap)
    GCl.Dec.value      = GCl.Dec + off_x* 146.48/(subp*GCl.cosmoPS*3600)/(1+GCl.mockobs.z_snap)
    GCl.massproxy_mask = mask
    
    """ Masks relic emission outside R200 of the given cluster """
    if image2mask is not None:
        masked_im, mask = Mask(image2mask, radius, (cen_x[0], cen_y[0]), setto=setto)
        return off_x, off_y, masked_im
    else:
        return off_x, off_y

    
def Mask(image, radius, center, setto=-10):
    """ masks in image outside the influence zone of a given cluster based on the mass """

    # compare maximum with local maximum in Rvir --> is there a more massive component?
    # find minimal contour within rvir that divides both clusters
    # use this contour to mask the outer regions
    # another massive cluster in the image
    mask      = cmask(image,radius, center=center)  #(1+GCl.z)
    masked_im = np.copy(image)
    masked_im[mask] = setto
#    plt.imshow(masked_im)
#    plt.show()
    
    return masked_im, mask




def weight_snap(snaplist_z, z, sigma_z=0.2, use_list=None, fake=False):
    """ Takes:
        snaplist_z : A list of usable snaps  
                 z : The redshift of the assigned cluster
                 
        Returns: An numpy array of the weight for each snapshot
    """
    
    # Draws the gaussian distribution function
    if fake:
        return np.ones_like(snaplist_z)
    
    vetted = gaussian_pseudo(np.asarray([z-zz for zz in snaplist_z]), 0, sigma_z)
    
    if use_list is None: # If not specified, make all snapshots available
        use_list = [True] * len(snaplist_z)
    
    # Weights snaps according to their gaussian weight
    weighted_snap = np.multiply(vetted, np.asarray(use_list))
    weighted_snap /= np.sum(weighted_snap)  
                      
    return weighted_snap


def assign_snaps(snaplist_z, boundaries_z, VCM, snaps_Nclusters, sigma_z=0.2, use_list=None, skycoverage=1, Vsimu=0.1,
                 fake=False, logmode='short', verbose=False):
    """ Takes:
        snaplist_z     : A list of usable snaps  
        boundaries_z   : The cental redshift of the volume
        VCM            : The comoving volume within the shell [Gpc**3]
        snaps_Nclusters: The number off individual clusters in each snapshot
        rand           : a random intitalizer for reproductive results
        
        Returns        : trials a list of lists that tells you about the exact number of representations for each cluster in each snap
    In the future 'sigma' should be given for the age as it is a more linear representation for the evolution of galaxy clusters in the universe
    
    Currently this is only done for one z value. The goal would be to create a list for an array of z values
    + initial randomization that then could be saved
    """
    
    clusterlist = []
    z_central = np.mean(boundaries_z)
    weighted_snap = weight_snap(snaplist_z, z_central, sigma_z=sigma_z, use_list=use_list)
    weighting = VCM * skycoverage * weighted_snap/Vsimu

    # Now do a poisson trial for each cluster in each snapshot
    #print(snaplist_z)
    #print(weighting * snaps_Nclusters[0])
    #return weighting * snaps_Nclusters[0]

    for kk, (weight, N_clusters) in enumerate(zip(weighting, snaps_Nclusters)):
         if not fake:
             trials = np.random.poisson(weight, N_clusters)
         else:
             trials = np.ones(N_clusters)
         
         # do while there are clusters in the trial; exspected are natural number array
         while np.sum(trials) > 0:
             positions          = np.where(trials > 0)
             trials[positions] -= 1
             clusterlist       += zip([kk]*len(positions[0]), positions[0])

    if verbose:
        if logmode == 'long':
            print('z_central %.4f' % z_central, 'weighted_snap', weighted_snap, 'p_counts',
                  ['%.1e' % (w*N) for w,N in zip(weighting,snaps_Nclusters)])
        elif logmode == 'short':
            print('z_central %.4f' % z_central, 'N_GCl', len(clusterlist))

    return clusterlist         # return all clusters within shell


def discovery_prop(relics, a=3, b=0.10):
    """ Takes a list of relics and takes returns an array of they weighted discovery probabilities according 
    to their shape parameter and the discovery function 
    
    Input: relics a list of CLusterBuster relics
           a is the normalization of the strength
           b is the offset
           """
    from scipy.stats import logistic
    
    probs = logistic.cdf([np.log10(b/relic.shape_advanced().value)*a for relic in relics])
    return probs


#   Take from [zstart to zend] take deltaV or delta Z --> Compute deltaV and mutliply with completeness sky  to get observed volume
#   Use observed volume and 
