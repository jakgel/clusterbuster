# not working, I really have to learn how to install python modules
# this set path question makes me ccraysi!
#import os
#import sys
#sys.path.insert(0, os.getcwd())
from astropy.io import fits  #import pyfits as fits
import warnings
import numpy as np
import matplotlib.pyplot as plt



#====== Hickjacking the header from another NVSS image  
#hduTemplate = fits.open('./../pyutil_my/Template-FITS/Casa4.3_Image.fits')
hduTemplate = fits.open('./../pyutil_my/Template-FITS/NVSS_Image.fits') 
#hdulist.info()
hduTemplateHead = hduTemplate[0]    
     
     
     
def FITS2numpy (fitsfile) :    
  
  hdulist = fits.open(fitsfile)

  hdu        = hdulist[0]
  hduHead    = hdu.header 

  image      = hdu.data.squeeze() #hdu.data[0,0,:,:]
  pixelsize  =  ( hduHead['CDELT1'], hduHead['CDELT2'] )
  center = [ [hduHead['CRVAL1'],hduHead['CRVAL2']] , [hduHead['CRPIX1'],hduHead['CRPIX2']] ]
  
  #k = 0
  #image = np.swapaxes(np.swapaxes(image, 0, k)[::-1], 0, k)



  return (image, center, pixelsize)  
     
# Writing FITS
def numpy2FITS ( array,  outfile, spixel, center=[-1,-1], pcenter = [-1,-1], oname='', nuobs=0, telescop='UnknownRadioTelescope') : #, header=hdutemplateHead.header

    ''' Creates an FITS file out of an 2D numpy map 
   
       input: 
      
       spixel (float) in arcsecs
       The created images ar eup to my knowledge not readable by CASA
    '''
  

    #print hduTemplateHead.header   
      
    #print header['NAXIS'] # = array.shape[0]
    #header['NAXIS1'] = array.shape[0]
    #header['NAXIS2'] = array.shape[1]
    #print header['NAXIS1'], header['NAXIS2'], header['NAXIS3'], header['NAXIS4']
     
    
    #NBINS = 500
    #histogram = plt.hist(array.flat, NBINS)
    #plt.show()
     
    newarray   = array[np.newaxis , np.newaxis, :, :].astype(np.float32) #array[np.newaxis , np.newaxis, :, :]

    hdu        = fits.PrimaryHDU(newarray)
#    hduHead    = hdu.header  # now add some modifications ...    
    hduHead    = hduTemplateHead.header   
     
    # Debugging delete certai entries
    keylist = ['PC01_01', 'PC02_01', 'PC03_01', 'PC04_01',
               'PC01_02', 'PC02_02', 'PC03_02', 'PC04_02',
               'PC01_03', 'PC02_03', 'PC03_03', 'PC04_03',
               'PC01_04', 'PC02_04', 'PC03_04', 'PC04_04']
#    for key in keylist:
#        try: 
#            del hduHead[key]
#        except:
#            print '[%s] not found, so we cannot delete it from the .fits header' % (key)

    #NAXIS   =                    2   #4                                                  
    hduHead['NAXIS1']  =   array.shape[0]                                                  
    hduHead['NAXIS2']  =   array.shape[1]                                                  
    #NAXIS3  =                    1                                                  
    #NAXIS4  =                    1 
 
    hduHead['EXTEND']  =          True #T  
    hduHead['CDELT1']  =  -1*float(spixel)/3600    #[deg]; if not float will set 3600   
    hduHead['CDELT2']  =   1*float(spixel)/3600    #[deg]; if not float will set 3600  
    hduHead['HISTORY'] = ''
    
    if center[0] != -1:
      hduHead['CRVAL1']  =   float(center[0])          #RA --SIN                                                                                                                           
      hduHead['CRVAL2']  =   float(center[1])          #Dec--SIN          
      hduHead['OBSRA']   =   float(center[0])                                               
      hduHead['OBSDEC']  =   float(center[1])      
      
    if pcenter[0] != -1:
      hduHead['CRPIX1']  =   float(pcenter[0]) 
      hduHead['CRPIX2']  =   float(pcenter[1])  
     
    #hduHead['BSCALE']  =   1.000000000000E+00 /PHYSICAL = PIXEL*BSCALE + BZERO                 
    #hduHead['BZERO']   =   0.000000000000E+00                                                  
    #hduHead['BMAJ']    =   2.772961722480E-03                                                  
    #hduHead['BMIN']    =   2.552167309655E-03                                                  
    #hduHead['BPA']     =  -4.343997192383E+01                                                  
    #hduHead['BTYPE']   = 'Intensity'                                                           
    hduHead['OBJECT']  = oname                                                     
                                                      
    #hduHead['BUNIT']   = 'Jy/beam '           /Brightness (pixel) unit                         
    #hduHead['EQUINOX'] =   2.000000000000E+03                                                  
    #hduHead['RADESYS'] = 'FK5     '                                                            
    #hduHead['LONPOLE'] =   1.800000000000E+02                                                  
    #hduHead['LATPOLE'] =   3.798332668224E+01                                                  
    #PC01_01 =   1.000000000000E+00                                                  
    #PC02_01 =   0.000000000000E+00                                                  
    #PC03_01 =   0.000000000000E+00                                                  
    #PC04_01 =   0.000000000000E+00                                                  
    #PC01_02 =   0.000000000000E+00                                                  
    #PC02_02 =   1.000000000000E+00                                                  
    #PC03_02 =   0.000000000000E+00                                                  
    #PC04_02 =   0.000000000000E+00                                                  
    #PC01_03 =   0.000000000000E+00                                                  
    #PC02_03 =   0.000000000000E+00                                                  
    #PC03_03 =   1.000000000000E+00                                                  
    #PC04_03 =   0.000000000000E+00                                                  
    #PC01_04 =   0.000000000000E+00                                                  
    #PC02_04 =   0.000000000000E+00                                                  
    #PC03_04 =   0.000000000000E+00                                                  
    #PC04_04 =   1.000000000000E+00                                                  
    #hduHead['CTYPE1']  = 'RA---SIN'                                                            
    #hduHead['CRVAL1']  =   2.057500144906E+02                                                                                                 
    #hduHead['CRPIX1']  =   3.201000000000E+03                                                  
    #hduHead['CUNIT1']  = 'deg     '                                                            
    #hduHead['CTYPE2']  = 'DEC--SIN'                                                            
    #hduHead['CRVAL2']  =   3.798332668224E+01                                                                                              
    #hduHead['CRPIX2']  =   3.201000000000E+03                                                  
    #hduHead['CUNIT2']  = 'deg     '                                                            
    #hduHead['CTYPE3']  = 'FREQ    '                                                            
    #hduHead['CRVAL3']  =   1.7E+09                                                  
    #hduHead['CDELT3']  =   5.0E+07                                                  
    #hduHead['CRPIX3']  =   1.000000000000E+00                                                  
                                                    
    #hduHead['CTYPE4']  = 'STOKES  '                                                            
    #hduHead['CRVAL4']  =   1.000000000000E+00                                                  
    #hduHead['CDELT4']  =   1.000000000000E+00                                                  
    #hduHead['CRPIX4']  =   1.000000000000E+00                                                  
    #hduHead['CUNIT4']  = '        '                                                            
    #hduHead['PV2_1']   =   0.000000000000E+00                                                  
    #hduHead['PV2_2']   =   0.000000000000E+00   
    
    if nuobs > 0:
        hduHead['CUNIT3']  = 'GHz      '  
        hduHead['RESTFRQ'] = 1.4E+09  #1.7E+09  /Rest Frequency (Hz)               #3.250000000000E+08                    
    #hduHead['SPECSYS'] = 'TOPOCENT'           /Spectral reference frame                        
    #hduHead['ALTRVAL'] =  -7.034932238758E+01 /Alternate frequency reference value             
    #hduHead['ALTRPIX'] =   1.000000000000E+00 /Alternate frequency reference pixel             
    #hduHead['VELREF']  =                  259 /1 LSR, 2 HEL, 3 OBS, +256 Radio                 
    #hduHead['COMMENT'] casacore non-standard usage: 4 LSD, 5 GEO, 6 SOU, 7 GAL                 
    #hduHead['TELESCOP']= 'Hoeft+2015 SIMULATION    '                                                            
    #hduHead['OBSERVER']= 'C. KONAR'                                                            
    #hduHead['DATE-OBS']= '2003-09-21T10:14:10.369263'                                          
    #hduHead['TIMESYS'] = 'TAI     '                                                                                                         
    #hduHead['OBSGEO-X']=   1.657069564687E+06                                                  
    #hduHead['OBSGEO-Y']=   5.797969000072E+06                                                  
    #hduHead['OBSGEO-Z']=   2.073044091645E+06                                                  
    #hduHead['DATE']    = '2015-02-25T18:27:09.044000' /Date FITS file was written              
    #hduHead['ORIGIN']  = 'CASA 4.4.85 (DEV r85)'                                          
    
    hduHead['OBSERVER'] = 'MockObs'
    hduHead['TELESCOP'] = telescop
    del hduHead['HISTORY'] # Just to much of false history ... lets delete it
    
    
    with warnings.catch_warnings(): # To catch the warning if an image is overwritten
       warnings.simplefilter('ignore')
       fits.writeto( filename=outfile, data=newarray, header = hduHead, clobber = True, checksum=False)
       
    return outfile

#====
import scipy.ndimage.morphology    as morp #morp.binary_dilation(input, structure=None, iterations=1, mask=None, output=None, border_value=0, origin=0, brute_force=False)
import pylab 


def numpy2mask ( array, cutoff, Ndil, outfile=False) : #, header=hdutemplateHead.header
  
 array2 = np.greater(array, (array*0.)+cutoff)
 array3 = morp.binary_dilation(array2, iterations=Ndil)
 
 # call numpy2FITS ( array,  outfile) with bool as data type
 newarray   = array3[np.newaxis , np.newaxis, :, :].astype(np.float32) #np.bool_) #array[np.newaxis , np.newaxis, :, :]

 hdu        = fits.PrimaryHDU(newarray)
 hduHead    = hdu.header  # now add some modifications ...
     
 hduHead    = hduTemplateHead.header   
                                                  
 hduHead['NAXIS1']  =   array.shape[0]                                                  
 hduHead['NAXIS2']  =   array.shape[1]                                                  
 
 hduHead['EXTEND']  =     True           #    T               
 if outfile:
  fits.writeto( filename=outfile, data=newarray, header = hduHead, clobber = True)
    
 return array3
 
 