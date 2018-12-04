#!/usr/bin/env python


"""
Created on 2015

@author: jakobg


 This script is about extracting (NVSS)-relics when giving it an properly formated .csv file of galaxy clusters. It expects for every galaxy cluster
 - one .fits file in the 'Images_NVSS' folder
 - one .region file in the 'Regions' folder which contouns at least one DS9 region called polygon with text describing the regions name (e.g. drection), relic type and spectral index
 [optionally] one .slist file in the 'Sources' folder for source subtraction before extracting relics
 output is an .csv file of all extracted relics
 
"""


from __future__ import division,print_function

import os
import math
import numpy  as np
import pandas as pd
import argparse

    
print( '###==== Step 0a: Executing self written .py subroutines ====###' )
import clusterbuster.surveyclasses    as cbclass
import clusterbuster.sourceextraction as relex
import clusterbuster.dbclasses        as cdb
import clusterbuster.iout.misc        as iom
import clusterbuster.iout.surveyplots as ioclass
import clusterbuster.maput            as maput
import clusterbuster.fitsut           as fitsut


from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from astropy.wcs import WCS
from reproject import reproject_interp
                
def updateClusters_missingRegions(ClList, AddList):
    """ This adds relic regions based on a .csv file to the cluster List 
    It is a workaround to have a fast update of the already existing datafiles with the soon to be standard file system,
    in this case  only for clusters without NVSS detections are updated via the identifiers.
    
    Later on this might be incorporated in the CB suite
    """

    regions = pd.read_csv(AddList, comment='#')
    """ FIX the classification by mapping"""    
    mapping = {'PHOENIX'          : 0, 
               'RELIC'            : 1,
               'SHOCKLET'         : 3,
               'AGN'              : -1,
               'AGN_relic'        : -1,  
               'AGN_RELIC/PHOENIX': 0,
               'HALO'             : -2,
               'SHOCKLET'         : 3}
    
    regions["CLASS_int"] = regions["FLAGS_CLASS"].map(mapping)
    regions["CLASS_int"] = regions["CLASS_int"].fillna(-1)
    regions.loc[regions['FLAGS_CLASS'] == 1 & regions['Counter'], 'FLAGS_CLASS'] = 2      # Double Relic
    regions.loc[regions['FLAGS_CLASS'] == 3 & regions['Counter'], 'FLAGS_CLASS'] = 2      # Double Relic
     

    for index, row in regions.iterrows():
        """ Just append the cluster if its status is not True and add the missing relic regions """
        for GCl in ClList:
            if GCl.status not in ['TRUE'] and GCl.name == row['Cluster']:  

                if int(row['CLASS_int']) == -2:
                    continue  # For now exclude halos, just because they are currently a seperate entry in the clusters file
                rtype = int(row['CLASS_int'])
                alpha = -1  #row['Alpha']
                alpha_err = 0  #row['Alpha_error']
                alphaFLAG = False
                candidate = (row['FLAG_conf'] == False)
                region = cbclass.RelicRegion(name=row['Identifier'], cnt=[], rtype=rtype, alpha=alpha, alpha_err=alpha_err,
                                             alphaFLAG=alphaFLAG, candidate=candidate)
                GCl.add_regions([region])
                
    return ClList


def survey_run(surveys, infolder='', outfoldertop='/data/ClusterBuster-Output/', plot=True):
    """ Extracts survey relics from an real world survey
    """
    
    for survey in surveys:
        print('###==== Step 0b: Initialize internal variables/objects for survey: %s ====###' % survey)
        smt = iom.SmartTiming()
        
        ClList = []
        Excluded = []
        subtract = ['slist', 'fits']  # im, fits, slist

        Jy_SI    = 1e-26    # W/Hz/m^2
        outfolder = '%s%s' % (outfoldertop, survey)
        topfolder = os.getcwd() # '/home/jakobg/lib/ClusterBuster/Relics_Surveys/'
        iom.check_mkdir(outfolder)  # create folder if necesairy

        print( '###==== Step 1: Load data and anaylise it   ====###' )
        # np.genfromtxt('Analysis_RORRS/ClusterRelics.csv'', delimiter=';')
        ClusterFile = infolder + 'ClusterList/ClusterAfterNuza2017_clusters.csv' 
        RegionFile  = infolder + 'ClusterList/ClusterAfterNuza2017_regions.csv' 
        
        
        Clusters = pd.read_csv(ClusterFile, comment='#', delimiter=',', quotechar='"')
        
        Clusters.where(Clusters.notnull(), 0)
        
        """ Part of development: rpelace nan values with values that can be handled by clustebruster """
        
        for strings in ['REF_LX', 'REF_M200', 'REF_M500', 'REF_F']:
            Clusters[strings] = Clusters[strings].replace(np.nan, '', regex=True)

        for values in ['M200', 'M500', 'LX_500_0.2-2.4']:
            Clusters[values] = Clusters[values].replace(np.nan, 0, regex=True)

        n = 0
        for index, CL in Clusters.iterrows():
            if CL['Cluster'] and CL['Cluster'] not in [o.name for o in ClList] and CL['Cluster'] not in ['']:
                """I did this to remove unfinished, but recent additions to the relic database"""
                #print(type(CL['Discovery']))
                #if math.isnan(CL['Discovery']):
                #    pass
                #elif int(CL['Discovery']) >= 2018:
                #    continue

                try:
                    """I did this to remove unfinished, but recent additions to the relic database"""
                    if CL['Cluster'] == '#':
                        continue
                    if int(CL['Discovery']) >= 2018:
                        continue
                except:
                    pass


                n += 1
#                if n > 5:
#                     continue
                Cl_name = CL['Cluster']
                status  = CL['FLAG_INCLUDED']
                

                print(CL)

                RA_host  = float(CL['RA'])    #float(CL[2])  
                Dec_host = float(CL['Dec'])
                
                diff = np.sqrt((float(CL['RA_Xmax'])-RA_host)**2+(float(CL['Dec_Xmax'])-Dec_host)**2)*3600
                if not math.isnan(diff) :           
                    print('Two different centre positions given for cluster %s. The offset was %.1f arcsec' % (Cl_name, diff) )
                    RA_host  = float(CL['RA_Xmax'])
                    Dec_host = float(CL['Dec_Xmax'])
  
                z    = float(CL['z'])
                M200 = float(CL['M200'])*1e14
                M500 = float(CL['M500'])*1e14
                Lx   = float(CL['LX_500_0.2-2.4'])
                flux_lit = float(CL['F_lit'])
                  
                halo = CL['FLAG_Halo'] 
              
                try:
                  ClassFlag = ('true' in CL['Type_Flag'].lower())
                except :
                  ClassFlag = False    
                  
                #create Class object
                GCl = cbclass.Galaxycluster(name=Cl_name, RA=RA_host, Dec=Dec_host, z=z, M200=M200, M500=M500, Lx=Lx, 
                                            Lx_lit=Lx, flux_lit=flux_lit, ClassFlag=ClassFlag, halo=halo, status=status)

                # add further references
                GCl.Lx      .ref = cdb.reference(CL['REF_LX'  ], rtype='text', page=None              , nr=None)
                GCl.M200    .ref = cdb.reference(CL['REF_M200'], rtype='text', page=CL['REFPAGE_M200'], nr=None)
                GCl.M500    .ref = cdb.reference(CL['REF_M500'], rtype='text', page=CL['REFPAGE_M500'], nr=None)
                GCl.flux_lit.ref = cdb.reference(CL['REF_F'   ], rtype='text', page=CL['REFPAGE_F']   , nr=None)


                #============= Load  survey (NVSS) image  =============#
                if GCl.status not in ['TRUE']: 
                    ClList.append(GCl)
                    continue
                fitsimage = infolder + 'Images_%s/%s-%s.fits' % (survey, survey, Cl_name)
                image, center, spixel = fitsut.fits2numpy(fitsimage)
                
                
                if survey == 'NVSS':
                        s_pixel   = [spixel[1]*GCl.cosmoPS*3600, spixel[1]*3600]
                        NVSSbeam  = [45., 45./s_pixel[1]]
                        NVSS_rms  = 4.5e-4     # in Jy/beam
                        NVSSlimit = 2*NVSS_rms
                        NVSSnu    = 1.4
                        telescope = 'VLA-D'
                        GCl.dinfo = cbclass.DetInfo(beam=[NVSSbeam[0], NVSSbeam[0], 0],
                                                     spixel=s_pixel[1], rms=NVSS_rms, 
                                                     limit=NVSSlimit, telescope=telescope, 
                                                     nucen=NVSSnu, center=center[0], pcenter=center[1])
                if survey == 'TGSS':
                        s_pixel   = [spixel[1]*GCl.cosmoPS*3600, spixel[1]*3600]
                        TGSSbeam  = [25., 25./s_pixel[1]]
                        TGSS_rms  = 3.0e-3     # in Jy/beam
                        TGSSlimit = 2*TGSS_rms
                        beamrec   = 1. if GCl.Dec > 19 else 1. / np.cos(np.radians(GCl.Dec - 19))
                        TGSSnu    = 0.1475
                        telescope = 'GMRT'
                        GCl.dinfo = cbclass.DetInfo(beam=[TGSSbeam[0]*beamrec, TGSSbeam[0], 0],
                                                     spixel=s_pixel[1], rms=TGSS_rms, 
                                                     limit=TGSSlimit, telescope=telescope, 
                                                     nucen=TGSSnu, center=center[0], pcenter=center[1])

                #============= Load relic search region  =============#
                # Make in np.image
                regfile     = infolder + 'Regions/RR_%s.reg' % (Cl_name) 
                GCl.regions = ioclass.readDS9relics(regfile, spixel, center[0], center[1])

                #============= Subtract Sources  =============#
                # in Sources folder
                #try load folder:
                    #img = bdsm.process_image(args.file+args.ft, thresh_isl=args.tIs, thresh_pix=args.tPi, mean_map = 'zero', beam = (0.0125,0.0125,0), rms_map = False, rms_value = 0.00045, thresh = 'hard') 
                #except:

                #pybdsm.catalog_type
                #--< create .fits image ut of that, which you subtract from your image ....
                smt(task='subtraction')
                model            = np.zeros((image.shape))
                model_conv       = np.zeros((image.shape))
                use_list, use_im = (False, False)
                if 'slist' in subtract:
                    slist = infolder + 'Sources/slist/%s.slist' % Cl_name
                    if os.path.isfile(slist):
                        scL = iom.read_para_list( slist )

                        for sc in scL:

                            if sc['shape'] == 'Gaussian':
                                g_size = [float(sc['majoraxis'])/s_pixel[1]*60, float(sc['minoraxis'])/s_pixel[1]*60]
                            else:
                                g_size = [GCl.dinfo.beam[0]/GCl.dinfo.spixel, GCl.dinfo.beam[1]/GCl.dinfo.spixel]
                            freq_factor = (GCl.dinfo.nucen/1.4)**(-0.7)
                            COOp = iom.CoordinateToPixel(iom.J2000ToCoordinate(sc['dir']), spixel, center[0], center[1])

                            #This is not good --> better create an unconcolved model and convolve it with the desired beam
                            model += maput.ImageGaussian_inv(model.shape, sc['flux']*1e-3*freq_factor, g_size, [COOp[0]-1.,COOp[1]-1], theta = sc['theta'], FWHM=True)  #*gaussian_area
                            #model_conv += maput.ImageGaussian_inv(model_conv.shape, sc['flux']*1e-3*freq_factor, g_size, [COOp[0]-1.,COOp[1]-1], theta = sc['theta'], FWHM=True)  #*gaussian_area
                        model_conv = model
                        use_list = True

                if 'fits' in subtract:
                    highres_image = infolder + 'Images_%s/%s-%s.fits' % ("FIRST", "FIRST", Cl_name)
                    if os.path.isfile(highres_image):
    
                        # regridd
       
                        # http://reproject.readthedocs.io/en/stable/  --> works on fits files
                        hdu1 = fits.open(fitsimage)[0]
                        image_2, center_2, spixel_2 = fitsut.fits2numpy(highres_image)
                        s_pixel_2    = [spixel_2[1]*GCl.cosmoPS*3600,spixel_2[1]*3600]  
                        fitsut.numpy2fits ( image_2 ,  infolder + 'Images_%s/%s-%s_test.fits' % ("FIRST", "FIRST", Cl_name), s_pixel_2[1], center_2[0], center_2[1]) 
                        hdu2 = fits.open(infolder + 'Images_%s/%s-%s_test.fits' % ("FIRST", "FIRST", Cl_name))[0]      
                        hdu2.data = hdu2.data.squeeze()
                        hdu2.data[np.isnan(hdu2.data)]      = 0.     # For contour masked  NVSS images I encountered the isue that some values where nan
                        hdu2.data[np.where(hdu2.data<6e-4)] = 0.
                                  
                        pad = 50
                        hdu2.data = np.lib.pad(hdu2.data, pad, maput.padwithtens)         
                        
                        FWHM2sigma     =  1/2.354 
                        FWHM_FIRST     =  5.4
                        FWHM_conv      = np.sqrt(GCl.dinfo.beam[0]**2-FWHM_FIRST**2)
                        gaussian_2D_kernel = Gaussian2DKernel(FWHM_conv/s_pixel_2[1]*FWHM2sigma)  
                        A_beam_old     =  1.133*((FWHM_FIRST/s_pixel_2[1])**2)  # FIRST-beam
                        A_beam         =  1.133*((GCl.dinfo.beam[0]/s_pixel_2[1])**2)
                        from copy import copy
                        hdu2_conv      =  copy(hdu2)                 
                        hdu2_conv.data = A_beam/A_beam_old*convolve(hdu2.data, gaussian_2D_kernel, normalize_kernel=True) 
                        for hdu in [hdu2, hdu2_conv]:
#                            hdu.data = np.expand_dims(hdu.data, axis=0)
#                            hdu.data = np.expand_dims(hdu.data, axis=0)
                            hdu.header['CRPIX1'] = hdu.header['CRPIX1'] + pad
                            hdu.header['CRPIX2'] = hdu.header['CRPIX2'] + pad

#
#                        from astropy.io import fits
#                        from astropy.utils.data import get_pkg_data_filename
#                        hdu1 = fits.open(get_pkg_data_filename('galactic_center/gc_2mass_k.fits'))[0]
#                        hdu2 = fits.open(get_pkg_data_filename('galactic_center/gc_msx_e.fits'))[0]
###                      
                        for hdu in [hdu1,hdu2_conv]:
#                            print(hdu.header)
                            try:
                                hdu.data = hdu1.data[0,0,:,:]
                            except:
                                print('Test ... data dimensions are matching')
                            hdu.header['NAXIS']      = 2

                            keylist = ['PC01_01', 'PC02_01', 'PC03_01', 'PC04_01',
                                       'PC01_02', 'PC02_02', 'PC03_02', 'PC04_02',
                                       'PC01_03', 'PC02_03', 'PC03_03', 'PC04_03',
                                       'PC01_04', 'PC02_04', 'PC03_04', 'PC04_04',
                                       'NAXIS3','NAXIS4',
                                       'CTYPE3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CROTA3',
                                       'CTYPE4','CRVAL4','CDELT4','CRPIX4','CUNIT4','CROTA4']
                            for key in keylist:     
                                try: 
                                    del hdu.header[key]
                                    print('[%s] removed from the .fits header' % (key))
                                except:
                                    print('[%s] not found, so we cannot delete it from the .fits header' % (key))
                                    
                            hdu.header['EPOCH']     = 3e3
                            hdu.header['EQUINOX']   = 3e3
                            print('====================') 
                            
                        fitsut.numpy2fits ( hdu2_conv.data.squeeze(),  infolder + 'Images_%s/%s-%s_test2.fits' % ("FIRST", "FIRST", Cl_name), s_pixel_2[1], center_2[0], [c+pad for c in center_2[1]])  


                        print( 'WCS(hdu1.header).wcs.naxis, WCS(hdu2.header).wcs.naxis', WCS(hdu1.header).wcs.naxis, WCS(hdu2_conv.header).wcs.naxis )
    
       
#                        array, footprint = reproject_interp(hdu2, hdu1.header) #hdu 2 image and systm, hdu1--> just the system
                        array, footprint = reproject_interp(hdu2_conv, hdu1.header) #hdu 2 image and systm, hdu1--> just the system 
#                        array, footprint = reproject_exact(hdu2_conv, hdu1.header) #does not solve the problem                      
                        print('_______', np.sum(array), footprint)
                       

                        print('_______________________________', array.shape, image.shape)  
                        array = array.squeeze() #could be removed

                        fitsut.map2fits(array , GCl.dinfo, infolder + 'Images_%s/%s-%s_test3.fits' % ("FIRST", "FIRST", Cl_name))    
                        squeezed = array.squeeze() #could be removed
                        squeezed[np.isnan(squeezed)]      = 0.
                        print( 'fits_subtraction: np.sum(squeezed):', np.sum(squeezed), ', np.sum(array.squeeze()):', np.sum(array.squeeze()) )
                        
                        model_conv  = squeezed # add up  OR replace!   
                        use_im      = True
#                    except:
#                      warnings.warn("No .fits specified for %s. Alternatively the format of the .fits file causes problems."  % (Cl_name))

                residuum  = image-model_conv

                """ Development: Only get the flux within the search region """       
                extreme_res =  True
                residuum =  maput.ContourMasking(residuum,[region.cnt[0] for region in GCl.regions])
                

                
                print( '%30s source subtraction;  list: %5r; image: %5r'   % (Cl_name , use_list, use_im) )
                GCl.maps_update(residuum  , 'Diffuse'    , infolder + '%s/Images_%s/diffuse/%s-%s.fits'            % (topfolder, survey, survey, Cl_name))
                if np.sum(model_conv) != 0 or extreme_res:      
                    GCl.maps_update(image     , 'Raw'       , infolder + '%s/Images_%s/raw/%s-%s_res.fits'         % (topfolder, survey, survey, Cl_name))
                    GCl.maps_update(model     , 'Modell'    , infolder + '%s/Images_%s/subtracted/%s-%s.fits'      % (topfolder, survey, survey, Cl_name))
                    GCl.maps_update(model_conv, 'Subtracted', infolder + '%s/Images_%s/subtracted/%s-%s_conv.fits' % (topfolder, survey, survey, Cl_name))
                smt()
                
                #============= impose relic.search  =============#
                for ii,region in enumerate(GCl.regions):
                  smt(task='RelicExtr')
                  relics = relex.RelicExtraction(residuum, z, GCl=GCl, dinfo=GCl.dinfo, rinfo = region, Imcenter=center, subtracted=model)[0:2]  # faintexcl=3.6
                  smt()                                
                  relics = sorted(relics, key=lambda x: x.flux, reverse=True)
                  
                  for relic in relics:
                      relic.alpha.value = region.alpha
                  
                  GCl.add_relics(relics)   

                # Add galaxy cluster to the list
                ClList.append(GCl)
          
            #============= Report why certain clusters are excluded  =============#
            else:
                RL_name = CL['Cluster']
                if CL['Identifier']:
                    RL_name += '_' + CL['Identifier']

                if CL['FLAG_INCLUDED'] in ['noMAP']:
                    string = RL_name + ' excluded because the corresponding region is not mapped by the survey.'
                else:
                    string = RL_name + ' excluded because of: ' + CL['FLAG_INCLUDED']

                Excluded.append(string)

            mf = open("%s/Excluded.dat" % outfolder, "w")
            for ex in Excluded:
                mf.write(ex + '\n')
                 
        ClList = sorted(ClList, key=iom.Object_natural_keys)

    
        print('#=====  Last Step: Output is produced ====#'); smt(task='output') 
        
        """ This is an intervening step: Update and ... the missing clusters, in the future this might done at the beginning at an first step """
        ClList = updateClusters_missingRegions(ClList, RegionFile)  # topfolder+RegionFile
        

        print('#=====  A: Pickle Objects ====#' )
        iom.pickleObject(ClList, outfolder+'pickled/', 'ClList')
        
        print('#=====  B: Create the Survey and pickle it ====#')
        cnt_levels = [9e-4, 1.8e-3, 3.6e-3, 7.2e-3, 1.44e-2]

        synonyms = [('1RXS J060313.4+421231', '1RXS J06+42'),
                    ('ACT-CLJ0102-4915', 'ACT-CLJ01-49'),
                    ('CIZA J0649.3+1801', 'CIZA J0649'),  # CIZA J0649+18
                    ('CIZA J0107.7+5408', 'CIZA J0107'),
                    ('CIZA J2242.8+5301', 'CIZA J2243'),  # CIZA J2243+53
                    ('MACS J0025-1222', 'MACS J0025'),
                    ('MACS J0717.5+3745', 'MACS J0717'),  # MCS J0717+37
                    ('MACS J1149.5+2223', 'MACS J1149'),  # J1149+22
                    ('MACS J1752.0+4440', 'MACS J1752'),  # MCS J1752+44
                    ('MACS J2243.3-0935', 'MACS J2243'),  # MCS J1752+44
                    ('MaxBCG 138.91895+25.19876', 'MaxBCG 138+25'),
                    ('MaxBCG 217.95869+13.53470', 'MaxBCG 217+13'),
                    ('PSZ1 G004.5-19.5', 'PSZ1 G004'),  # PSZ1 G004-19
                    ('PSZ1 G096.89+24.17', 'PSZ1 G097'),  # PSZ1 G097+24
                    ('PSZ1 G108.18-11.53', 'PSZ1 G108'),  # PSZ1 G108-12
                    ('PLCK G200.9-28.2', 'PLCK G200'),  # PSZ1 G108-12
                    ('PLCK G287.0+32.9', 'PLCK G287'),  # PLCK G287+33
                    ('RXC J0225.1-2928', 'RXC J0225'),
                    ('RXC J1053.7+5452', 'RXC J1054'),  # RXC J1054+55
                    ('RXC J1053.7+5452 ', 'RXC J1053'),  # RXC J1054+55
                    ('RXC J1234.2+0947', 'RXC J1234'),  # RXC J1054+55
                    ('RXC J1314.4-2515', 'RXC J1314'),  # RXC J1314-25
                    ('ZwCl 0008+5215', 'ZwCl 0008'),  # ZwCl 0008+52
                    ('ZwCl 1447.2+2619', 'ZwCl 1447'),  # ZwCl 0008+52
                    ('ZwCl 2341+0000', 'ZwCl 2341'),  # ZwCl 2341+00
                    ('[KMA2007] 217.95869+13.53470', 'KMA2007'),  # ZwCl 2341+00
                    ]

#        synonyms_lit = [('2017A&A...597A..15D', '2017A+A_deGasperin+Intema+'),
#                        ('2017arXiv170801718K', '2017arXiv_Kale+Wik+')]
        for GCl in ClList:

            for syn in synonyms:
                if GCl.name.replace('_', ' ') == syn[0] or GCl.name == syn[0]:
                    print( '! Name replacement:', syn)
                    GCl.name = syn[1]
            GCl.name = GCl.name.replace('_', ' ')
            GCl.updateInformation()

        norm = cdb.norm('R200', Nexp=1)
        Histo = cdb.Histogram2D(nbins=(64, 46), fromto=[[0, 2.*np.pi], [0, 1.5]], norm=norm)     # angle_projected(rad), D_proj(R200)
        Survey = cbclass.Survey(ClList, survey, cnt_levels=cnt_levels, synonyms=synonyms, dinfo=GCl.dinfo, hist_main=Histo, surshort=survey)  # 'NVSS' should be replaced with a real survey class
        Survey.emi_max = 2e-2
        Survey.scatterkwargs = {"alpha": 0.7, "fmt": "o", "markersize": 10}
        Survey.histkwargs = {"alpha": 0.4}
        iom.pickleObject(Survey, outfolder+'/pickled/', 'Survey')

        for GCl in Survey.GCls:
            print(GCl.name, GCl.status)

    smt(forced=True)
    return True


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Extracts survey relics from an real world survey')
    parser.add_argument('-surveys', dest='surveys', nargs='+', action='store', default=['NVSS'],  type=str, help='Survey names that are to be used')
    parser.add_argument('-infolder', dest='infolder', action='store', default='', type=str, help='filepath for data arrays to be loaded (folder of NVSS paths etc.)')
    parser.add_argument('-outputfolder', dest='outputfolder', action='store', default='/data/ClusterBuster-Output/', type=str, help='filepath for data arrays to store')
    args = parser.parse_args()

    survey_run(surveys=args.surveys, infolder=args.infolder, outfoldertop=args.outputfolder, plot=False)
