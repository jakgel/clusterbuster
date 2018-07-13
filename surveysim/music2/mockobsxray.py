#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 17:26:48 2017

@author: jakobg
"""

#!/usr/bin/env python
# Version JG --> 
# Remark: In order to simulate a radio intoferometric  observation it will be neccesairy to at least for some examples create visibilities and apply a NVSS deconvolution routine to them. 
# It is important to specify this so that we can be sure that to effect is important (or not!)

# Known BUGS:




from __future__ import division,print_function
import numpy                as np
import surveysim.music2.loadsnap   as LSconv
#import CreateRadioCubes     as crrc

import scipy.ndimage               as ndi 
import clusterbuster.mathut        as maut
import clusterbuster.iout.misc     as iom
import clusterbuster.constants     as myu
import clusterbuster.surveyclasses as cbclass

#from clusterbuster.XRayAnalysis    import *


#from pyutil.load_snap       import *

def Run_MockObs_XRay( gcl, savefolder, log=False, PCA=[], saveFITS = True, headerc='Mpc', verbose = True):
    ''' Creates maps of projected quantities in the galaxie cluster, like
          X-ray 
          Mass 
          Velocity field 
    '''
    
    
    
    mockobs    = gcl.mockobs
    z          = gcl.z.value
    
    smt = iom.SmartTiming(rate=5e4) #, logf=outf+'smt.log'  #;  print '###==== Step2a:  Loading configuration files ====###'  
  

    strSn = ('%s/SHOCKS_%05i/cluster.%05i.snap.%03i.shocks') % (gcl.mockobs.xrayfolder, gcl.mockobs.clid, gcl.mockobs.clid, gcl.mockobs.snap)   
    if verbose: print('Loading snapshot:',strSn)
    snap    = LSconv.Loadsnap(strSn,headerc=gcl.mockobs.headerc)  
    

  
    hSize       =  mockobs.hsize                       # kpc     
    HSize       =  mockobs.hsize*snap.head['hubble']   # kpc/h   multiplied with 0.7  ... reduced Hubble constant used in MUSIC-2

    rot         = maut.Kep3D_seb(Omega=mockobs.theta, i=mockobs.phi,omega=mockobs.psi)
    posrot      = np.dot( snap.pos, rot )
    velrot      = np.dot( snap.vel, rot )

    smt(task='PCA_[sub]')
    xSize     = 0.8*hSize 
    XSize     = 0.8*HSize 
    XSize_z  = XSize*(1+z)   # This comes fr
        
#    if xPixel*gcl.cosmoPS < 10:    # If the cluster is very nearby, use an less resolved  X-Ray/Mass image for the analysis
#      xPixel    = 10/gcl.cosmoPS   
      
    '''DEVELOPMENT: uses an very well resolved X-ray image '''
    xPixel    = 10/gcl.cosmoPS  #Size of pixel in arcseconds, choosen so that the resolution is 10 kpc
    
    nbins_x   =  int(2 * xSize / (gcl.cosmoPS*xPixel) )
#            x_pixel   = [xPixel*gcl.cosmoPS, xPixel]
    Xbeam     = [4*xPixel, 4]
    A_Xbeam   = 1.133*Xbeam[1]*Xbeam[1]
    mapinfo   = cbclass.DetInfo(beam= [0, 0, 0],  
                                spixel= xPixel, rms=0, 
                                limit = 0, telescope='MUSIC-2 simulation', 
                                nucen = 0, center=[0,0], pcenter=[nbins_x/2,nbins_x/2])
    
    
    if verbose: print('  ###==== Step 3a:  Binning cluster data cube  (X-Ray) ====###')

#    iL  =  np.where( snap.rho[:] >  -1 )[0] #snap.xray[:] >  -1
    
    # iL  =  np.where( ( -hThick < posrot[:,2] ) & ( posrot[:,2] < hThick )  &  ( snap.radi[:] > -1e60 ) )


    # bremsstrahlung, compare https://www.mrao.cam.ac.uk/~kjbg1/lectures/lect3.pdf

    fac   = LSconv.conversion_fact_gadget_rho_to_nb( snap.head )*LSconv.conversion_fact_ne_per_nb()
    brems = np.power(snap.rho*fac,2)*np.power(snap.u,1/2) # mimics the total emissivity of bremsstrahlung; neclects metallicity; neclects that MUSIC-2 doesn't have cooling
    iL  =   np.where( ( (np.log10(snap.rho*fac)+5) -(np.log10(snap.u)-2.2)*(3.8/2.8) ) < 0)[0] #snap.rho < 5e-4 with inverse color scalling           
   

    '''  BREMSTRAHLUNG & MASS
         (snap.rho < 3e-1) & (snap.u > 3e4) to remove the galaxies 
         (snap.rho < 5e-4) for a nice contrast image
         
         np.where( (snap.rho < 3e-1) & (snap.u > 3e4) )[0] for galaxy exclution
         
         np.where( (snap.rho > 3e-1) & (snap.u < 3e4) )[0] for galaxies only
         
                   
       For division galaxies and ICM, ...
       log10(rho)/log10(u)
       -5/2.2
       3.8 2.8
       -1.2/5
       
        np.where( ( (np.log10(snap.rho)+5) -(np.log10(snap.u)-2.2)*(3.8/2.8) ) < 0)[0] for galaxies only
         
    '''
    H_bems, xedges, yedges = np.histogram2d ( -posrot[iL,0], -posrot[iL,1], weights=    brems[iL] , range=[[-XSize_z,XSize_z], [-XSize_z,XSize_z]], bins=nbins_x )  #
    H_mass, xedges, yedges = np.histogram2d ( -posrot[iL,0], -posrot[iL,1],                         range=[[-XSize_z,XSize_z], [-XSize_z,XSize_z]], bins=nbins_x ) #weights=np.ones( (iL.shape[0])), 
    H_bems_conv = A_Xbeam*ndi.gaussian_filter(H_bems, (myu.FWHM2sigma*Xbeam[1],myu.FWHM2sigma*Xbeam[1]))  ## gaussian convolution, factor 0.5 steems from the fact that beam with is two times the gaussian standard-deviation
    H_mass_conv = A_Xbeam*ndi.gaussian_filter(H_mass, (myu.FWHM2sigma*Xbeam[1],myu.FWHM2sigma*Xbeam[1]))  ## gaussian convolution, factor 0.5 steems from the fact that beam with is two times the gaussian standard-deviation           
  
    
         
    ''' VELOCITY FIELD
    Here we derive the quantities needed to plot a mass weighted velocity field
    '''                                   
    H_gas_x, xedges, yedges = np.histogram2d ( -posrot[iL,0], -posrot[iL,1], weights=velrot[iL,0], range=[[-XSize_z,XSize_z], [-XSize_z,XSize_z]], bins=nbins_x )
    H_gas_x_conv            = A_Xbeam*ndi.gaussian_filter(H_gas_x, (myu.FWHM2sigma*Xbeam[1],myu.FWHM2sigma*Xbeam[1]))                                            
  
    # Get the real average velocity, Subtract the average central velocity
    H_gas_x                = H_gas_x     /H_mass      - H_gas_x_conv[int(nbins_x/2),int(nbins_x/2)]/H_mass_conv[int(nbins_x/2),int(nbins_x/2)]
    H_gas_x_conv           = H_gas_x_conv/H_mass_conv - H_gas_x_conv[int(nbins_x/2),int(nbins_x/2)]/H_mass_conv[int(nbins_x/2),int(nbins_x/2)]
    
    H_gas_y, xedges, yedges = np.histogram2d ( -posrot[iL,0], -posrot[iL,1], weights=velrot[iL,1] , range=[[-XSize_z,XSize_z], [-XSize_z,XSize_z]], bins=nbins_x ) 
    H_gas_y_conv            = A_Xbeam*ndi.gaussian_filter(H_gas_y, (myu.FWHM2sigma*Xbeam[1],myu.FWHM2sigma*Xbeam[1]))                                                                      
      
    H_gas_y                = H_gas_y     /H_mass      -  H_gas_y_conv[int(nbins_x/2),int(nbins_x/2)]/H_mass_conv[int(nbins_x/2),int(nbins_x/2)]
    H_gas_y_conv           = H_gas_y_conv/H_mass_conv -  H_gas_y_conv[int(nbins_x/2),int(nbins_x/2)]/H_mass_conv[int(nbins_x/2),int(nbins_x/2)]
                                       
#    H_gas_speed             = np.sqrt(H_gas_x     **2 + H_gas_y     **2)            
    H_gas_speed_conv        = np.sqrt(H_gas_x_conv**2 + H_gas_y_conv**2)            
    H_gas_speed_conv[np.where(H_mass_conv < 1e2)] = 0
    
    H_gas_angle             = np.arctan( H_gas_y     /H_gas_x)        * 180 / np.pi
    H_gas_angle_conv        = np.arctan( H_gas_y_conv/H_gas_x_conv)   * 180 / np.pi
                                          
                                             
                                             
    if saveFITS:
       iom.check_mkdir(savefolder) 
       fitstypes  = ['Brems', 'MassSpeed', 'MassAngle','Mass']    
       fitsarray = [np.log10(H_bems_conv).clip(min=-6.5)+6.5,H_gas_speed_conv,H_gas_angle_conv,np.log10(H_mass_conv)]  # np.clip(np.log10(H_bems_conv), 0, -9)
       
       for IM,fitstype in zip(fitsarray,fitstypes):
           fitsname = '%s/maps/z%04.f/%s-%s.fits' % (savefolder, gcl.mockobs.z_snap*1000, gcl.name, fitstype)
           if verbose: print('Gonna save', fitsname)
           gcl.maps_update(IM,  fitstype, fitsname, dinfo=mapinfo)


       H_x = H_gas_angle#indarr[0,:] + 1  #add sdsdsdsds
       H_y = H_gas_angle#indarr[1,:] + 1  #dssdsds
                      
       fitstypes = ['Ra', 'Dec', 'dx','dy' , ]    
       fitsarray = [H_x,H_y,H_gas_x_conv,H_gas_y_conv]  # np.clip(np.log10(H_bems_conv), 0, -9)
       
       
       print(gcl.mapdic['Brems'])
       
       for IM,fitstype in zip(fitsarray,fitstypes):
           fitsname = '%s/maps/z%04.f/%s-%s.fits' % (savefolder, gcl.mockobs.z_snap*1000, gcl.name, fitstype)
           if verbose: print('Gonna save', fitsname)
           gcl.maps_update(IM,  fitstype, fitsname, dinfo=mapinfo)
    return 
