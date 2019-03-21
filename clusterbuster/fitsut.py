""" A two way road for a conversion from numpy arrays to -fits files with the correct format
    Currently the header is initialized by a template.fits which should be changed in the future!
    I could also provide the funtionality to provide CASA (standard radioastronomical software) readable images
"""


from __future__ import division,print_function
from astropy.io import fits

import os
import warnings
import clusterbuster
import numpy as np
import scipy.ndimage.morphology as morp


#====== Hickjacking the header from another NVSS image  
#hduTemplate = fits.open('./../pyutil_my/Template-FITS/Casa4.3_Image.fits')
path = os.path.dirname(clusterbuster.__file__)
print(path)
hduTemplate = fits.open(path+'/Template-FITS/NVSS_Image.fits') 
#hdulist.info()
hduTemplateHead = hduTemplate[0]    
     

def fits2numpy (fitsfile):
  
  hdulist = fits.open(fitsfile)

  hdu = hdulist[0]
  hduHead = hdu.header

  image = hdu.data.squeeze() #hdu.data[0,0,:,:]
  pixelsize  =  ( hduHead['CDELT1'], hduHead['CDELT2'] )
  center = [ [hduHead['CRVAL1'],hduHead['CRVAL2']] , [hduHead['CRPIX1'],hduHead['CRPIX2']] ]

  return (image, center, pixelsize)  

def map2fits (array,dinfo,outfile):
    """
    Maps numpy info and detinfo directly to outfile
    This function uses some of the ObjectClass.py data, so it should in principle become part it.
    Detection info should be changed to mapinfo + detinfo (detection threshold etc.) 
    
    dinfo: ClusterBuster 
    
    """
    numpy2fits(array, outfile, dinfo.spixel, center=dinfo.center, pcenter=dinfo.pcenter, nuobs=dinfo.nucen,
               beam=dinfo.beam, telescope=dinfo.telescope)
    

def numpy2fits(array, outfile, spixel, center=None, pcenter = None, oname='', nuobs=1.4, beam=[45/3600,45/3600,0],telescope='UnknownRadioTelescope'):

    """ Creates an FITS file out of an 2D numpy map 
   
       input: 
      
       spixel (float) in arcsecs
       The created images are up to my knowledge not readable by CASA
       centre  = [Ra,Dec]
       pcenter = [p1,p2]
       
       Possible improvement: directly supply it as a function of a map class
       
    """

    newarray = array[np.newaxis , np.newaxis, :, :].astype(np.float32) #array[np.newaxis , np.newaxis, :, :]

    """=== This is a big designe desicion: Do I want to stick with two our four dimension. 
    Four dimension seem to be needed for the astrply imager and cas. Two dimensions are needed for the rotation packgae """
    hdu = fits.PrimaryHDU(newarray)
    hduHead = hdu.header  # now add some modifications ...
#    hduHead    = hduTemplateHead.header   

    hduHead['CTYPE1'] = 'RA---SIN'
    hduHead['CTYPE2'] = 'DEC--SIN'
    
    #NAXIS   =                    2   #4
    hduHead['NAXIS1'] = array.shape[0]
    hduHead['NAXIS2'] = array.shape[1]
    #NAXIS3  =                    1
    #NAXIS4  =                    1
 
    hduHead['EXTEND'] =  True
    hduHead['CDELT1'] = -1*float(spixel)/3600    # [deg]; if not float will set 3600
    hduHead['CDELT2'] =  1*float(spixel)/3600    # [deg]; if not float will set 3600
    hduHead['HISTORY'] = ''
    
    if center is not None:
      hduHead['CRVAL1'] = float(center[0])          # RA --SIN
      hduHead['CRVAL2'] = float(center[1])          # Dec--SIN
      hduHead['OBSRA'] = float(center[0])
      hduHead['OBSDEC'] = float(center[1])
      
    if pcenter is not None:
      hduHead['CRPIX1'] = float(pcenter[0])
      hduHead['CRPIX2'] = float(pcenter[1])
      
    hduHead['BUNIT'] = 'JY/BEAM '
    hduHead['BMAJ'] = beam[0]/3600  # .beam.set_major(GCl.dinfo.beam[0] * u.arcsecond)
    hduHead['BMIN'] = beam[1]/3600  #  f.beam.set_minor(GCl.dinfo.beam[1] * u.arcsecond)
    hduHead['BPA']  = beam[2]       # f.beam.set_angle(GCl.dinfo.beam[2])  # degrees

    hduHead['OBJECT'] = oname
    hduHead['EPOCH'] = (2000, 'Celestial coordinate equinox')

    if nuobs > 0:
        hduHead['RESTFRQ'] = (nuobs*1e9, 'in Hz')

    hduHead['OBSERVER'] = 'MockObs'
    hduHead['TELESCOP'] = telescope
    del hduHead['HISTORY'] # Just to much of false history ... lets delete it
    
    
    with warnings.catch_warnings(): # To catch the warning if an image is overwritten
       warnings.simplefilter('ignore')
       fits.writeto( filename=outfile, data=newarray, header = hduHead, overwrite = True, checksum=False)
       
    return outfile

def numpy2mask(array, cutoff, Ndil, outfile=False):
  
    array2 = np.greater(array, (array*0.)+cutoff)
    array3 = morp.binary_dilation(array2, iterations=Ndil)

    newarray = array3[np.newaxis , np.newaxis, :, :].astype(np.float32) #np.bool_) #array[np.newaxis , np.newaxis, :, :]

    hdu     = fits.PrimaryHDU(newarray)
    hduHead = hdu.header  # now add some modifications ...
    hduHead = hduTemplateHead.header

    hduHead['NAXIS1'] = array.shape[0]
    hduHead['NAXIS2'] = array.shape[1]
    hduHead['EXTEND'] = True

    if outfile:
        fits.writeto( filename=outfile, data=newarray, header = hduHead, clobber = True)

    return array3
 


def sparse_array_to_fits(GCls, outfolder, maptype = "Diffuse", source_type="relics"):

    for GCl in GCls:
        kpc_to_arcsec = 1/(GCl.cosmoPS)

        bins = ((np.asarray(range(int(GCl.dinfo.pcenter[0]*2+1))) - GCl.dinfo.pcenter[0])+0.25) * GCl.dinfo.spixel  # np.linspace(-500,500,widht=1)
        array_full = np.zeros((int(GCl.dinfo.pcenter[0]*2), int(GCl.dinfo.pcenter[1]*2)))

        if source_type == "relics":
            expression = lambda x: x.relics
        elif source_type == "compacts":
            expression = lambda x: x.compacts

        if maptype == "Temperature":
            sparseWeight = lambda x: x.wT
        elif maptype == "Density":
            sparseWeight = lambda x: x.wRho_up
        elif maptype == "Mach":
            sparseWeight = lambda x: x.wMach
        else:
            sparseWeight = lambda x: x.sparseW


        for relic in expression(GCl):

            x = np.sin(relic.sparseA) * relic.sparseD * kpc_to_arcsec
            y = np.cos(relic.sparseA) * relic.sparseD * kpc_to_arcsec

            array, x, y = np.histogram2d(x, y, bins=[bins, bins], weights=sparseWeight(relic))

            array_full += array

        mapname = maptype
        print('%s/maps/%s/%s.fits' % (outfolder, maptype.lower(), GCl.name))
        fitsname = '%s/maps/%s/%s.fits' % (outfolder, maptype.lower(), GCl.name)
        GCl.maps_update(array_full, mapname, fitsname)