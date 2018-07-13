#!/usr/bin/env python
# Version JG --> 
# Remark: In order to simulate an radio intoferometric  observation it will be neccesairy to at least for some examples create visibilities and apply a NVSS deconvolution routine to them. 
# It is important to specify this so that we can be sure that to effect is important (or not!)

# Known BUGS:

from __future__ import division,print_function

print( '###==== Step 1: Executing .py subroutines====###' )

import numpy        as np
import copy

import surveysim.music2.loadsnap     as loadsnap

import clusterbuster.surveyclasses    as cbclass
import clusterbuster.mathut           as mathut
import clusterbuster.iout.misc        as iom
import clusterbuster.maput            as maput
import clusterbuster.sourceextraction as relex
import clusterbuster.constants        as myu

#import clusterbuster.load_snap                 as LSconv
#from clusterbuster.XRayAnalysis    import *



#from pyutil.load_snap       import *

def SPH_binning(snap, posrot, dinf, iL, immask=1, HSize_z=1000, nbins=100, weights=None, convolve=True):
    
    H, xedges, yedges = np.histogram2d( -posrot[iL,0], -posrot[iL,1], weights=weights(snap), range=[[-HSize_z,HSize_z], [-HSize_z,HSize_z]], bins=nbins) #norm=LogNorm(), ,  cmin=1e-3   
    
    if convolve:
       H =  dinf.convolve_map(H*immask)                     
                                   
    return H


    

def Run_MockObs( bulked, GClrealisations, locations, steps = [1,2,3,4,5,6,7], CASAmock=False, XRay=False, saveFITS = False, 
                writeClusters = False, savewodetect=False, log=False, tPi=2, iii=0, sideEffects=False):
    ''' Runs a mock observation
        sideEffects: put True if you want the input galaxy clsuter to be changed, False if you want only a copy to be influenced '''
    (snap, Rmodel, emptySurvey) = bulked    
    efflist     = Rmodel.effList


    savefolder = emptySurvey.outfolder
    iom.check_mkdir(savefolder) 

    #1. Python: Create convolved perfect simulational output image
    #2. PyBDSM/Python: Create clean mask on that with some dilation
    #3. Casa/Python: Create constant rms and beam-corrected image (clean)
    #4. Casa/Python: Apply this clean mask with immath on constant-rms image
    #5. PyBDSM/python: Use pybdsm with detection_image= 'masked_constant rms_image'
    #6. Python. Create masked .fits mock imagename
    #7  Python. Extract radio relics

    #Variant B: Clean mask and .fits --> Source parameters; like variant A from step 4 on
    if CASAmock:
        import drivecasa as drica
        casa = drica.Casapy()
    
   
    smt = iom.SmartTiming(rate=5e4) #, logf=outf+'smt.log'  #;  print( '###==== Step2a:  Loading configuration files ====###'  )
  
  
    #  Units, conversion factors, and input variables
    fac_rho      =  loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb() #electrons/cm^-3
    fac_T        =  loadsnap.conversion_fact_gadget_U_to_keV(  snap.head ) # in [keV]                   
    fac_T2       =  loadsnap.conversion_fact_gadget_U_to_keV(  snap.head )/8.61732814974056e-08  # to K          
                           
    ''' determines if you want to change the galaxy cluster or not '''
    if sideEffects:  GClrealisations_used =                GClrealisations       
    else          :  GClrealisations_used =  copy.deepcopy(GClrealisations)           
    for jj, gcl in enumerate(GClrealisations_used):
        mockobs    = gcl.mockobs
        z          = gcl.z.value
        hsize      = mockobs.hsize
        dinf       = gcl.dinfo  #copy.deepcopy(gcl.dinfo) # A local realisation of dinfo. This is just because some parameters of dinfo could change, because of adaptive pixelsize etc.
        fac_x      = loadsnap.comH_to_phys( snap.head, z )
        
        if log: print( '  Realisation #%i-%i for %i efficiencies. In total %i realisations' %(iii+1, jj+1, len(efflist), len(GClrealisations) ) )
        #  Load variables and setting survey parameters
        for ii,eff in enumerate(efflist):
       
            if ii==0 and log:   
                print( 'Begin with efficiency of %5.3e.' %(eff) )
          
    
            #  Units, conversion factors, and input variables
            radiounit   = myu.radiounit_A*eff # erg/s/Hz    --- Unit of particle luminousity in .radio snaps
            rot         = mathut.Kep3D_seb(Omega=mockobs.theta, i=mockobs.phi,omega=mockobs.psi)
            posrot      = np.dot( snap.pos, rot ) * fac_x
            #velrot      = np.dot( snap.vel, rot )  Taken out, as long as we don't need to plot the velocity vetors

            
            smt(task='Bin_radio_[sub]'); #print( '###==== Step 3b:  Binning cluster data cube  (radio) ====###'
            # Parameters implied
            # See Nuza+ 2012 Equ (1)
            relativistics  = (1+z)   # Term A: Bandwidth quenching  # Term B: Decreased Photon Emission Rate & Term C decreased photon energy both are included in gcl.cosmoDL
    #        s_radio     =  radiounit / Jy    / (4*np.pi* gcl.cosmoDL**2)      * relativistics   # radiounit*s_radio   is Jy/particle        #Umrechnung
            s_radio_SI  =  radiounit / myu.Jy_SI / (4*np.pi*(gcl.cosmoDL*1e-2)**2)* relativistics   # radiounit*s_radioSI is Jy/particle        #Umrechnung
            nbins       =  int(2 *hsize / (gcl.cosmoPS*dinf.spixel) )  
            if nbins>mockobs.binmax and ii==0:
                
                binsold       = nbins
                spixelold     = dinf.spixel   
                dinf.spixel   = dinf.spixel    * np.power(float(nbins)/float(mockobs.binmax), 0.5)
                mockobs.hsize = mockobs.hsize  * np.power(float(nbins)/float(mockobs.binmax),-0.5)
                dinf.update_Abeam()
                nbins         = mockobs.binmax
                if log: print( 'At z=%.3f with an pixel size of %.1f arcsec, the number of pixels per image is %i^2. Due to that the pixelscale was increased to %.1f arcsec and the binned boxsize decreased to %i kpc.' %  (z, spixelold, binsold, dinf.spixel, mockobs.hsize) )
                hsize         = mockobs.hsize
            dinf.pcenter = [nbins/2,nbins/2]
            
            ''' Initial: np.where( (snap.rdow < 5e-1/(loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb())) & 
                                   (snap.udow > 5e4) )[0] 
                         --> yielded a lot of compact emission associated to AGN
                         
            
                         np.where( (snap.rdow < 1e-1/(loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb())) 
                                  & (snap.udow > 1e5) 
                                  & (snap.mach < 8  ) )[0] 
                         --> still a lot of compact shock objects, but a better relic ratio
            
                        np.where( (snap.rdow < 5e-2/(loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb())) 
                              & (snap.udow > 3e5) 
                              & (snap.mach < 5  ) )[0] 
                        --> provided nearly a good result
            
            
                       CUT 00: Seems to be to harsch as it excludes a large portion of high mach components in radio relics.
                       np.where( (snap.rdow < 1e-2/(loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb())) 
                              & (snap.udow > 1e6) 
                              & (snap.mach < 5  ) )[0] 
                       
                       CUT 01: Less harsch, not working on r, ...
                           
                       np.where( (snap.rho  < 1e-2/(loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb())) 
                              & (snap.u    > 1e6) 
                              & (snap.mach < 9  ) )[0] 
                       
                       CUT 02: based on the seperation between galaxies and the ICM
                       
                       For division galaxies and ICM, ...
                       rho/T
                       -5/2.2
                       3.8 2.8
                       -1.2/5
                       
                       CUT_03: Prsented to my working group
                       (snap.u*fac_T>0.6 ) &
                       (( (np.log10(snap.rho*fac_rho)+6) + (np.log10(snap.u*fac_T)-0)*(3.8/2.8) ) > 0) &
                       (np.log10(snap.rho*fac_rho) > -5.3) & 
                       (snap.mach   < 12)       
            
                better to find a condition based on rho/u, because they describe the particle properties (i.e. galaxies)
            '''
        
            iL     =  np.where(         
#                            
                            ( (8.9+  3.3 - np.log(snap.u*fac_T2) - 0.65*np.log10(snap.rho*fac_rho)) < 0) &
                            ( (8.9- 11.0 - np.log(snap.u*fac_T2) - 3.5*np.log10(snap.rho*fac_rho))  < 0) &
                            ( (8.9+  6.9 - np.log(snap.u*fac_T2) + 0.5*np.log10(snap.rho*fac_rho))  < 0) &
                            (snap.mach   < 10)  &
                            (np.sqrt(snap.pos[:,0]**2+snap.pos[:,1]**2+snap.pos[:,2]**2)*fac_x < 1.5*gcl.R200() ) 
                           )[0]
 
            if hasattr(snap,'radiPre'):
               if log: print('Run_MockObs:: Ratio of PREs to total emission', ((np.sum(snap.radiPre[iL]))/(np.sum(snap.radi[iL])+np.sum(snap.radiPre[iL]))) )
               snap.radi += snap.radiPre
        
            H1, xedges, yedges = np.histogram2d( -posrot[iL,0], -posrot[iL,1], weights=s_radio_SI*snap.radi[iL], range=[[-hsize,hsize], [-hsize,hsize]], bins=nbins) #norm=LogNorm(), ,  cmin=1e-3   
            ''' DEVELOPMENT BEGIN
                trying to manage a simple subtraction of compact sources'''   
    
            ''' Differnce of gaussians method 
            
            We do this iteratively three times to also remove those particles that where masked by other 
            bright particles before
            
            This methos is defines by
            
            thresh: A threshold for masking
            scale_1: Smaller scale in kpc
            scale_2: Larger  scale in kpc
            
            
            '''
            thresh  = 0.7
            scale_1 = 75  #100
            scale_2 = 300 #450
            
            DoG1_filter        = copy.deepcopy(dinf)
            DoG1_filter.beam   = [scale_1/gcl.cosmoPS,scale_1/gcl.cosmoPS,0]        
            DoG1_filter.update_Abeam()
            
            DoG2_filter        = copy.deepcopy(dinf)
            DoG2_filter.beam   = [scale_2/gcl.cosmoPS,scale_2/gcl.cosmoPS,0]        
            DoG2_filter.update_Abeam()
            
            DoG_mask           = np.ones_like(H1)
            for ii in range(3):
                H2_first           = DoG1_filter.convolve_map(H1*DoG_mask)  ## gaussian convolution
                H2_sec             = DoG2_filter.convolve_map(H1*DoG_mask)  ## gaussian convolution
    #            with np.errstate(divide='ignore', invalid='ignore'):
                DoG_rel            = np.divide(np.abs(H2_sec-H2_first)+1e-20,H2_sec+1e-20)
                DoG_mask[np.where(DoG_rel < thresh)] = 0          #0.03                                                   
            H2_first           = DoG1_filter.convolve_map(H1)  ## gaussian convolution
            H2_sec             = DoG2_filter.convolve_map(H1)  ## gaussian convolution                                
                         
                         
            H2                 = dinf.convolve_map(H1*DoG_mask)
#            print('____ Masked/Unmasked flux (mJy):  %6.3f %6.3f' % (np.sum(H2)/dinf.Abeam[0]*1000,s_radio_SI*np.sum(snap.radi[iL])*1000))
                 

            ''' DEVELOPMENT END'''
            
            smt(task='WriteDilMask_[sub]'); #print( '###==== -- 4b:  Writing dilated .mask  ====###'
    #        mask       =  maput.numpy2mask (H2, dinf.limit, Ndil) # outfile = outfile.replace('.fits','') + '_mask.fits', 
            
            #img = bdsm.process_image(filename= outfile+'_simple_conv.fits', thresh_isl=args.tIs, thresh_pix=args.tPi, mean_map = 'zero', beam = (0.0125,0.0125,0), rms_map = False, rms_value = 0.00045, thresh = 'hard') 
            #img.export_image(outfile= outfile+'_simple_conv.ismk.fits'    , img_type='island_mask', img_format='fits', mask_dilation=5, clobber=True) 
            #img.export_image(outfile= outfile+'_simple_conv.ismk.mask'    , img_type='island_mask', img_format='casa', mask_dilation=5, clobber=True) 
    
            
            if CASAmock:
                ''' Removed implementation
                Any synthetic observations would need to be run in python within this shell. So we have to somehow put it as an module into python '''
    	        
            else:
                if log: print( '  ###====          - Using the simple convolved image ====###' )
                IM0 = H2  #(fits.open(simpleconv))[0].data 
                
            smt(task='CreateMask_[sub]'); #print( '###==== Step 6:  Create masked .fits mock image ====###'
            IM1  = np.squeeze(IM0)  #!!! Here unmasked! ... np.multiply(np.squeeze(IM0), mask) 
            '''Outdated saveFITS, please update and put to end of procedure '''
            if saveFITS and CASAmock:
              maput.numpy2FITS (IM1,  'sim.vla.d.masked.fits', dinf.spixel) 
              
            smt(task='RelicExtr_[sub]'); #print( '###==== Step 7: Extract radio relics ====###'
    
            relics = relex.RelicExtraction(IM1, s_radio_SI, z, GCl=gcl, dinfo=dinf, eff=eff, rinfo=cbclass.RelicRegion('',[],rtype=1)) #, faintexcl=0.4, Mach=Hmach, Dens=Hdens, 
    
            smt(task='RelicHandling_[sub]')
            relics          = sorted(relics, key=lambda x: x.flux, reverse=True)
            gcl.add_relics(relics) 
            

    
    #        print(('... ... ...', eff, np.sum(snap.radi), np.sum(snap.radi[iL]), np.sum(H2))
            if savewodetect or len(relics) > 0: 
                
            
                if eff == efflist[0]:
                    if log: print( '  ++The brightest relic found has a flux density of %f mJy' % (relics[0].flux) )  #Could producese errors, once there is no relict in the list                       
                    iom.check_mkdir(savefolder + '/maps/z%04.f'  % (gcl.mockobs.z_snap*1000))
                    #if saveFITS: gcl.Image = IM1 # The masked image is pushed into the galaxy cluster class, just for one relic in the hole list!  
               
               
                    ''' I couldnt come up with something better to take the inverse '''
                    mask = np.ones(snap.rho.shape,dtype=bool)
               
               
                    '''DEVELOPMENT ... if pixel in highres dominates its surroundings...(50 percent in 5*5 grid) ...add?'''
               
                    mask[iL]=0
                    Subtracted, xedges, yedges =  np.histogram2d( -posrot[mask,0], -posrot[mask,1], weights=s_radio_SI*snap.radi[mask], range=[[-hsize,hsize], [-hsize,hsize]], bins=nbins) #norm=LogNorm(), ,  cmin=1e-3   
                    Subtracted                 +=  H1*(1-DoG_mask)
               
                    SubConv                    =  dinf.convolve_map(Subtracted)
 
                ''' New and cool part to derive additional relic information like average and  mach number and alpha, 
                We also get the emission weighted density, as this works only on the bright parts it is fine to work with the subset of particles, yeah!
                '''

                Hmach,Hrho,Htemp,Harea,Halpha,Hmag,Hpre       = None,None,None,None,None,None,None
#                Hrho_down,Hrho_up,Htemp_down,Htemp_up  = None,None,None,None,
                alpha_help            =  (snap.mach[iL]**2+1)/(snap.mach[iL]**2-1)
            
                Hmach     = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.mach[iL])
                Halpha    = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*alpha_help)
                Hrho      = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.rho[iL] *fac_rho)
                Htemp     = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.u[iL]   *fac_T)  
                Harea     = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.area[iL])  
                Hmag      = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.B[iL])  
                Hpre      = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radiPre[iL])                
#                Hrho_up   = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.rup[iL] *fac_rho)
#                Hrho_down = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.rdow[iL]*fac_rho)
#                Htemp_up  = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.uup[iL] *fac_T)  
#                Htemp_down= SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.udow[iL]*fac_T)  
#                Htemp_rat = SPH_binning(snap, posrot, dinf, iL, immask=H2_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: x.udow[iL]/x.uup[iL])  
#                Hrho_rat  = SPH_binning(snap, posrot, dinf, iL, immask=H2_mask,HSize_z=hsize, nbins=nbins, weights=lambda x: x.rdow[iL]/x.rup[iL])  

                allflux = np.asarray([])
                for relic in relics:    
                    relic.wMach       =  Hmach[relic.pmask]     
                    relic.wRho        =  Hrho [relic.pmask]  
                    relic.wT          =  Htemp[relic.pmask] 
                    relic.wArea       =  Harea[relic.pmask] 
                    relic.wAlpha      =  Halpha[relic.pmask]      
                    relic.wB          =  Hmag[relic.pmask]
                    relic.wPre        =  Hpre[relic.pmask]
#                    relic.wRho_up     =  Hrho_up[relic.pmask]  
#                    relic.wRho_down   =  Hrho_down[relic.pmask]  
#                    relic.wT_up       =  Htemp_up[relic.pmask]  
#                    relic.wT_down     =  Htemp_down[relic.pmask]  
                    
                    relic.wDoG_rel  =  DoG_rel[relic.pmask]  
                    
                    allflux = np.concatenate((relic.sparseW,allflux), axis=0)           
                    relic.averages_quantities()
                    
                '''Save maps'''
                allflux = allflux.flatten()    
                if saveFITS:
                    ''' Here the maps are already masked with the detection region '''
                    
                    smt(task='WriteFits_[writes,sub]')
                    if log: print( '###==== Step 4:  Preparing FITS file & folders ====###' )
                    
                    parlist = (savefolder, gcl.mockobs.z_snap*1000, gcl.name, Rmodel.id)
                    gcl.maps_update(H1                                            , 'Raw'      , '%s/maps/z%04.f/%s-%04i_native.fits'          % parlist)
                    gcl.maps_update(IM1                                           , 'Diffuse'  , '%s/maps/z%04.f/%s-%04i.fits'                 % parlist)
                    gcl.maps_update(Subtracted                                    , 'CompModell','%s/maps/z%04.f/%s-%04i_compact.fits'         % parlist)
                    gcl.maps_update(SubConv                                       , 'Subtrated', '%s/maps/z%04.f/%s-%04i_compactObserved.fits' % parlist)
                    if(len(relics) > 0):
                        gcl.maps_update(gcl.Mask_Map(Hmach,normalize=allflux,eff=eff) , 'Mach'     , '%s/maps/z%04.f/%s-%04i_mach.fits'            % parlist)
                        gcl.maps_update(gcl.Mask_Map(Hrho ,normalize=allflux,eff=eff) , 'Rho'      , '%s/maps/z%04.f/%s-%04i_rho.fits'             % parlist)
                        gcl.maps_update(gcl.Mask_Map(Htemp,normalize=allflux,eff=eff) , 'Temp'     , '%s/maps/z%04.f/%s-%04i_temp.fits'            % parlist)
                        gcl.maps_update(gcl.Mask_Map(Hmag ,normalize=allflux,eff=eff) , 'B'        , '%s/maps/z%04.f/%s-%04i_B.fits'               % parlist)
                        gcl.maps_update(gcl.Mask_Map(Hpre,normalize=allflux,eff=eff)  , 'PreRatio' , '%s/maps/z%04.f/%s-%04i_prerat.fits'          % parlist)
                    gcl.maps_update(DoG_rel                                       , 'DoG_rel'  , '%s/maps/z%04.f/%s-%04i_DoG_rel.fits'         % parlist)
                    gcl.maps_update(DoG_mask                                      , 'DoG_mask' , '%s/maps/z%04.f/%s-%04i_DoG_mask.fits'        % parlist)


                ''' PhD feature --> plot the DoF images in a subplot
##                import matplotlib.pyplot as plt   
##                with np.errstate(divide='ignore', invalid='ignore'):
##                    DoG_rel            = np.divide(np.abs(H2_sec-H2_first)+1e-20,H2_sec+1e-20)
##                pixR200 =        gcl.R200()/(gcl.cosmoPS*dinf.spixel)     
##                bou     =  gcl.R200()*1.5   # pixR200*1
##                f, ((ax1, ax2, ax3)) = plt.subplots(1, 3, figsize=(30,12)) #, sharey='row', sharex='col', sharey='row'
##                ax1.imshow( np.power((np.abs(H2_first-H2_sec)), 1/1 )     , extent=(-hsize, hsize, -hsize, hsize)) #-bou+cen
##                im2 = ax2.imshow( DoG_rel                                 , extent=(-hsize, hsize, -hsize, hsize), vmin=0.2, vmax=1 ) #-bou+cen           
##                XX,YY = np.meshgrid(xedges[0:-1]+0.5*gcl.cosmoPS*dinf.spixel,yedges[0:-1][::-1]+0.5*gcl.cosmoPS*dinf.spixel) #yedges[-1:0]
##                ax2.contour(XX, YY, DoG_mask,  colors='r', levels=[0.5])
##                ax3.imshow( np.power(dinf.convolve_map(H1*DoG_mask), 1/1 ), extent=(-hsize, hsize, -hsize, hsize)  ) #-bou+cen
##                ax1.set_xlim(-bou, bou)
##                ax1.set_ylim(-bou, bou)
##                ax2.set_xlim(-bou, bou)
##                ax2.set_ylim(-bou, bou) 
##                ax3.set_xlim(-bou, bou)
##                ax3.set_ylim(-bou, bou)
##                
##                ax1.set_title('DoG')
##                ax2.set_title('DoG/LowResImage + mask (contours)')
##                ax3.set_title('Filtered NVSS')
##                
##                print('CreateMokObs',  pixR200, gcl.R200(),dinf.spixel, gcl.cosmoPS)
##                circle1 = plt.Circle((0, 0), gcl.R200(), fill=False, color='w', ls='-')
##                circle2 = plt.Circle((0, 0), gcl.R200(), fill=False, color='w', ls='-')      
##                circle3 = plt.Circle((0, 0), gcl.R200(), fill=False, color='w', ls='-')   
##
##                ax1.add_artist(circle1)
##                ax2.add_artist(circle2)
##                ax3.add_artist(circle3)
##                
##                cax2 = f.add_axes([0.42, 0.12, 0.2, 0.03]) 
##                cb2  = f.colorbar(im2, format='%.2f', ticks=[0.0, 0.25, 0.5, 0.75, 1.0], cax = cax2, orientation="horizontal")  #label='average Mach', 
##                
##                plt.savefig('%s/%s-%04i_joined.png'        % (savefolder, gcl.name, Rmodel.id)) #dpi=400
##                plt.savefig('%s/%s-%04i_joined.pdf'        % (savefolder, gcl.name, Rmodel.id)) #dpi=400          
#                '''     '''
                
                gcl.add_relics(relics) 
                PhD feature end '''
            
            else:
                if log: print( "  --No relic detected at efficiency of %5.3e or below for realisation #%i-%i. Stop search." %(efflist[ii-1], iii+1, jj+1) )
                break 

    if writeClusters:
        '''This is here because some outputs get lost in a multiprocessing heavy input/output queue process'''
        for gcl in GClrealisations_used:
            filename =  'GCl-%05i'  % (gcl.mockobs.id)
            iom.pickleObject( (gcl, Rmodel), savefolder + '/pickled/', filename, append = False) 


    if log: print( 'Finished with all efficency values' )
    return (True, smt, GClrealisations_used,  Rmodel)