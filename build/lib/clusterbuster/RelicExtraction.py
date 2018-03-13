''' Script for blob detection originaly inspired by http://www.learnopencv.com/blob-detection-using-opencv-python-c/ at 13th of Juny 2015
    Now it returns clusterBuster relics and hsitogramms
'''


# Standard imports
import numpy as np
import cv2
#import pyfits as fits
#import scipy.ndimage.morphology    as morp #morp.binary_dilation(input, structure=None, iterations=1, mask=None, output=None, border_value=0, origin=0, brute_force=False)
from   scipy import ndimage                # used to compute the centroid
from   scipy import sparse


#print '###==== Step 1: Loading .py subroutines====###'
import pyutil_my.ClusterBuster_pythonClass as CBclass
import pyutil_my.NPimageutil as pyNPi

def RelicExtraction(image, sradio, z, contourMask = [], Imcenter=False, GCl=False, dinfo=False, rinfo= False, faintexcl=False, subtracted=False, term=False, eff=1, cnt_minlength=5, HistoH=0):
  
    ''' Extracts Relics from an numpy array
    HistoH: -1 of not derived, 0 if derived, 1 if derived in the fly
    
    See http://photutils.readthedocs.io/en/stable/api/photutils.detect_sources.html
    for possible alternative of implementation
    '''
    
    
    if not GCl:
       GCl      = CBclass.Galaxycluster('NoName', 0, 0, z)
    if not dinfo: #If no detection info specified, create a default one
       dinfo    = CBclass.DetInfo()
    if not rinfo: #If no relic region is specified, create a default one
       rinfo    = CBclass.RelicRegion(name='', cnt=contourMask, rtype=-1, alpha=1.2)
    if not Imcenter:
      Imcenter = [ [0,0], [image.shape[0]/2,image.shape[1]/2] ]  # [ [Astronomical Coords], [Pixelcoords] ]
        
  
    if term: print '###==== Sub-Step X.1: Executing RelicExtraction(2d-numpy array) ====###' 
    #=== Define Image parameters
    relics  = []


    ## Way I of extracting relics: Relic detection by simple thresholding
    (contours, hierarchy)  = pyNPi.GetThreshContours(image, dinfo.limit, rinfo.cnt)   
    
    
    #img = removingcontours(img)   
    for jj,cnt in enumerate(contours): 
      if term: print '%6i out of %6i: len(cnt)=%6i' % (jj+1, len(contours), len(cnt) )
      M = cv2.moments(cnt)
      if cnt.shape[0] < cnt_minlength or M['m00'] == 0 : # cnt >= 3 needed for moment analysis & cnt >= 5 for ellipse fitting
          #print 'Contour of', len(cnt), 'ergo less than', cnt_minlength, 'points identified. It will not be considered further.'
          continue
      else:
        #cx = int(M['m10']/M['m00'])
        #cy = int(M['m01']/M['m00'])
        
        # Minimal enclosing circle
        (x,y),radius = cv2.minEnclosingCircle(cnt)
        
        # mask image, so that every value outside of the contour is set to zero
        region =  pyNPi.ContourMasking(image,[cnt])

        # compute flux and luminosity weighted coordinaetes/distance to cluster center
        flux             = 1e3*np.sum(region)/dinfo.Abeam[0] #in mJy
        GCl.flux_vol     = 1e3*np.sum(image )/dinfo.Abeam[0] #in mJy
        #====
        if faintexcl and flux < faintexcl:
            continue;
	
        try:
            subregion = pyNPi.ContourMasking(subtracted,[cnt])
            flux_ps   = 1e3*np.sum(subregion)/dinfo.Abeam[0]
        except:
            flux_ps   = 0

        centroid         = ndimage.measurements.center_of_mass(region) 
        
        # compute moment of inertia and gaussian fit
        cov, theta_iner  = pyNPi.inertial_axis(region)[2:4]

        
        Dec         = Imcenter[0][1] + (centroid[0]-Imcenter[1][1])*dinfo.spixel/3600
        RA          = Imcenter[0][0] - (centroid[1]-Imcenter[1][0])*dinfo.spixel/3600/np.cos(Dec*np.pi/180)        #Please mind the additional factor due to the projection!
 
 
        theta_elong = (np.degrees(theta_iner)+180)%180       
        LAS         = 2*(radius)*(dinfo.spixel/60)
        area        = (dinfo.spixel/60)**2*cv2.contourArea(cnt)  
        
        
        
        region[np.isnan(region)] = 0
        rel             = np.where(region>0)

        if HistoH >= 0:   
          ''' Creates the necessary masks and array for the FLUX histogram '''
          ''' The same thing with Mach number and ... weighted histogram would be equally important '''
          RA_arr          = Imcenter[0][0] - ( rel[1] - Imcenter[1][0])*dinfo.spixel/3600/np.cos(Dec*np.pi/180)
          Dec_arr         = Imcenter[0][1] + ( rel[0] - Imcenter[1][1])*dinfo.spixel/3600
          distance        = np.sqrt( np.add( np.power(np.abs(GCl.RA.value-RA_arr)*np.cos(Dec_arr*np.pi/180),2), np.power(np.abs(Dec_arr-GCl.Dec.value),2) ) ) * GCl.cosmoPS *3600      #2D-array in kpc, np.linalg.norm( [GCl.RA-RA, Dec-GCl.Dec] ) * GCl.cosmoPS *3600               
          angle_hist      = np.mod( np.arctan2(  -(GCl.RA.value-RA_arr)*np.cos(Dec_arr*np.pi/180) , Dec_arr-GCl.Dec.value ) + 1./2.*np.pi, 2.*np.pi)   #2D-array in kpc, np.linalg.norm( [GCl.RA-RA, Dec-GCl.Dec] ) * GCl.cosmoPS *3600
          weights         = region[rel]
                     
          Dproj_pix       = np.average(distance, weights=weights)  # in kpc 
          sparseD         = distance   #sparseD.csr_matrix(histo)  
          sparseA         = angle_hist
          sparseW         = weights
        else:
            Dproj_pix     = 0
            sparseD       = None
            sparseA       = None
            sparseW       = None 

#            d['mynewkey'] = 'mynewvalue'
            
        relic       = CBclass.Relic(rinfo, dinfo, RA, Dec, LAS, area, GCl=GCl, alpha=rinfo.alpha, alphaFLAG=rinfo.alphaFLAG, theta_elong=theta_elong, 
                                    cov=cov, cnt=cnt, F=flux, F_ps=flux_ps, Dproj_pix = Dproj_pix, eff=eff, 
                                    sparseD=sparseD, sparseW=sparseW, sparseA=sparseA, pmask=rel)
        relics.append( relic  )
        
        
        ''' Development: NOT functional. This would be used to create the Polar histogramms on the fly and 
        delete all contentional information afterwards to save memory '''
        if HistoH >= 1:
             relics.createHisto()
             relic.deletesparsedata()
    
    return relics
