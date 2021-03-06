#!/usr/bin/env python

from __future__ import division,print_function
import copy
import numpy as np
import surveysim.music2.loadsnap as loadsnap
import clusterbuster.surveyclasses    as cbclass
import clusterbuster.mathut           as mathut
import clusterbuster.iout.misc        as iom
import clusterbuster.maput            as maput
import clusterbuster.sourceextraction as relex
import clusterbuster.constants        as myu

def SPH_binning(snap, posrot, dinf, iL, immask=1, HSize_z=1000, nbins=100, weights=None, convolve=True):
    
    H, xedges, yedges = np.histogram2d( -posrot[iL,0], -posrot[iL,1], weights=weights(snap), range=[[-HSize_z,HSize_z], [-HSize_z,HSize_z]], bins=nbins) #norm=LogNorm(), ,  cmin=1e-3   
    
    if convolve:
       H = dinf.convolve_map(H*immask)
                                   
    return H


def Run_MockObs(bulked, GClrealisations, CASAmock=False, saveFITS=False, writeClusters=False,
                savewodetect=False, log=False, side_effects=False, filter_sp_phase=False, extract_subtracted=True):
    """ Runs a mock observation
        side_effects: put   True if you want the input galaxy cluster to be changed,
                            False if you want only a copy to be influenced """
    (snap, Rmodel, emptySurvey) = bulked
    savefolder = emptySurvey.outfolder
    iom.check_mkdir(savefolder)

    #Variant B: Clean mask and .fits --> Source parameters; like variant A from step 4 on
    if CASAmock:
        import drivecasa as drica
        casa = drica.Casapy()

    smt = iom.SmartTiming(rate=5e4)  # logf=outf+'smt.log'  #;  print( '###==== Step2a:  Loading configuration files ====###'  )

    #  Units, conversion factors, and input variables
    fac_rho = loadsnap.conversion_fact_gadget_rho_to_nb(snap.head)*loadsnap.conversion_fact_ne_per_nb()  #electrons/cm^-3
    fac_T   = loadsnap.conversion_fact_gadget_U_to_keV( snap.head)  # in [keV]
    fac_T2  = loadsnap.conversion_fact_gadget_U_to_keV( snap.head) / 8.61732814974056e-08  # to K
                           
    """ determines if you want to change the galaxy cluster or not """
    if side_effects:
        GClrealisations_used = GClrealisations
    else:
        GClrealisations_used = copy.deepcopy(GClrealisations)

    for jj, gcl in enumerate(GClrealisations_used):
        #  Load variables and setting survey parameters
        mockobs = gcl.mockobs
        z       = gcl.z.value
        hsize   = mockobs.hsize
        dinf    = gcl.dinfo    # Some parameters of dinfo could change, because of adaptive pixelsize etc.
        fac_x   = loadsnap.comH_to_phys(snap.head, z)
        eff = Rmodel.effList[0]

        #  Units, conversion factors, and input variables
        radiounit = myu.radiounit_A*eff # erg/s/Hz    --- Unit of particle luminousity in .radio snaps
        rot       = mathut.Kep3D_seb(Omega=mockobs.theta, i=mockobs.phi, omega=mockobs.psi)
        posrot    = np.dot(snap.pos, rot) * fac_x
        #velrot   = np.dot(snap.vel, rot)  Taken out, as long as we don't need to plot the velocity vetors

        smt(task='Bin_radio_[sub]'); #print( '###==== Step 3b:  Binning cluster data cube  (radio) ====###'
        # Parameters implied
        # See Nuza+ 2012 Equ (1)
        relativistics = (1+z)
        s_radio_SI  = radiounit / myu.Jy_SI / (4*np.pi*(gcl.cosmoDL*1e-2)**2) * relativistics   # radiounit*s_radioSI is Jy/particle        #Umrechnung
        nbins       = int(2 *hsize / (gcl.cosmoPS*dinf.spixel))
        if nbins > mockobs.binmax:
            binsold       = nbins
            spixelold     = dinf.spixel
            dinf.spixel   = dinf.spixel   * np.power(float(nbins)/float(mockobs.binmax), 0.5)
            mockobs.hsize = mockobs.hsize * np.power(float(nbins)/float(mockobs.binmax), -0.5)
            dinf.update_Abeam()
            nbins = mockobs.binmax
            if log:
                print('At z=%.3f with an pixel size of %.1f arcsec, the number of pixels per image is %i^2. Due to that the pixelscale was increased to %.1f arcsec and the binned boxsize decreased to %i kpc.' %  (z, spixelold, binsold, dinf.spixel, mockobs.hsize) )
            hsize = mockobs.hsize
        dinf.pcenter = [nbins/2, nbins/2]

        if filter_sp_phase:
            """ Filteres the cooled particles that no longer belong to the hot-ICM"""
            iL = np.where((np.sqrt(snap.pos[:,0]**2+snap.pos[:,1]**2+snap.pos[:,2]**2)*fac_x < 2.0*gcl.R200()) &
                            ((8.9 +  3.3 - np.log(snap.u*fac_T2) - 0.65*np.log10(snap.rho*fac_rho)) < 0) &
                            ((8.9 - 11.0 - np.log(snap.u*fac_T2) - 3.50*np.log10(snap.rho*fac_rho)) < 0) &
                            ((8.9 +  6.9 - np.log(snap.u*fac_T2) + 0.50*np.log10(snap.rho*fac_rho)) < 0) &
                            (snap.mach < 10)
                           )[0]
        else:
            iL = np.where(np.sqrt(snap.pos[:,0]**2+snap.pos[:,1]**2+snap.pos[:,2]**2)*fac_x < 2.0*gcl.R200())[0]

        if hasattr(snap, 'radiPre'):
            if log:
                print('Run_MockObs:: Ratio of PREs to total emission',
                      (np.sum(snap.radiPre[iL]))/(np.sum(snap.radi[iL])+np.sum(snap.radiPre[iL])))

        H1, xedges, yedges = np.histogram2d(-posrot[iL,0], -posrot[iL,1], weights=s_radio_SI*snap.radi[iL],
                                            range=[[-hsize,hsize], [-hsize,hsize]], bins=nbins)
        """ Difference of gaussians method - accomplishing a simple subtraction of compact sources"
        
        We do this iteratively three times to also remove those particles that where shadowed by other 
        bright particles before
        
        This method is defines by
        
        thresh: A threshold for masking
        scale_1: Smaller scale in kpc
        scale_2: Larger  scale in kpc        
        """
        thresh = 0.75
        scale_1 = 20
        scale_2 = 60

        DoG1_filter      = copy.deepcopy(dinf)
        DoG1_filter.beam = [scale_1/gcl.cosmoPS, scale_1/gcl.cosmoPS, 0]
        DoG1_filter.update_Abeam()

        DoG2_filter      = copy.deepcopy(dinf)
        DoG2_filter.beam = [scale_2/gcl.cosmoPS, scale_2/gcl.cosmoPS, 0]
        DoG2_filter.update_Abeam()

        DoG_mask = np.ones_like(H1)
        for no_use in range(2):
            convolved_sigma1 = DoG1_filter.convolve_map(H1*DoG_mask)  ## gaussian convolution
            convolved_sigma2 = DoG2_filter.convolve_map(H1*DoG_mask)  ## gaussian convolution
            DoG_rel = np.divide(np.abs(convolved_sigma2-convolved_sigma1)+1e-20, convolved_sigma2+1e-20)
            DoG_mask[np.where(DoG_rel < thresh)] = 0.2*DoG_mask[np.where(DoG_rel < thresh)]
        #convolved_sigma1 = DoG1_filter.convolve_map(H1)  ## gaussian convolution
        #convolved_sigma2 = DoG2_filter.convolve_map(H1)  ## gaussian convolution

        H2 = dinf.convolve_map(H1*DoG_mask)
#            print('____ Masked/Unmasked flux (mJy):  %6.3f %6.3f' % (np.sum(H2)/dinf.Abeam[0]*1000,s_radio_SI*np.sum(snap.radi[iL])*1000))

        smt(task='WriteDilMask_[sub]'); #print( '###==== -- 4b:  Writing dilated .mask  ====###'
#        mask       =  maput.numpy2mask (H2, dinf.limit, Ndil) # outfile = outfile.replace('.fits','') + '_mask.fits',

        #img = bdsm.process_image(filename= outfile+'_simple_conv.fits', thresh_isl=args.tIs, thresh_pix=args.tPi, mean_map = 'zero', beam = (0.0125,0.0125,0), rms_map = False, rms_value = 0.00045, thresh = 'hard')
        #img.export_image(outfile= outfile+'_simple_conv.ismk.fits'    , img_type='island_mask', img_format='fits', mask_dilation=5, clobber=True)
        #img.export_image(outfile= outfile+'_simple_conv.ismk.mask'    , img_type='island_mask', img_format='casa', mask_dilation=5, clobber=True)

        if CASAmock:
            """ Removed implementation; original idea:
            #1. Python: Create convolved perfect simulational output image
            #2. PyBDSM/Python: Create clean mask on that with some dilation
            #3. Casa/Python: Create constant rms and beam-corrected image (clean)
            #4. Casa/Python: Apply this clean mask with immath on constant-rms image
            #5. PyBDSM/python: Use pybdsm with detection_image= 'masked_constant rms_image'
            #6. Python. Create masked .fits mock imagename
            #7  Python. Extract radio relics """
        else:
            if log: print('###====          - Using the simple convolved image ====###')
            IM0 = H2  #(fits.open(simpleconv))[0].data

        smt(task='CreateMask_[sub]'); #print( '###==== Step 6:  Create masked .fits mock image ====###'
        IM1 = np.squeeze(IM0)  #!!! Here unmasked! ... np.multiply(np.squeeze(IM0), mask)
        """Outdated saveFITS, please update and put to end of procedure """
        if saveFITS and CASAmock:
          maput.numpy2FITS(IM1, 'sim.vla.d.masked.fits', dinf.spixel)

        smt(task='RelicExtr_[sub]');

        relics = relex.RelicExtraction(IM1, z, GCl=gcl, dinfo=dinf, rinfo=cbclass.RelicRegion('', [], rtype=1)) #, faintexcl=0.4, Mach=Hmach, Dens=Hdens,

        smt(task='RelicHandling_[sub]')
        relics = sorted(relics, key=lambda x: x.flux, reverse=True)
        gcl.add_relics(relics)
        if savewodetect or len(relics) > 0:

            if log: print('  ++The brightest relic found has a flux density of %f mJy' % (relics[0].flux))  #Could producese errors, once there is no relict in the list
            iom.check_mkdir(savefolder + '/maps/z%04.f' % (gcl.mockobs.z_snap*1000))

            """ Part to derive additional relic information like average and  mach number and alpha. We also get the 
            emission weighted density, as this works only on the bright parts it is fine to work with the subset of 
            particles
            """
            alpha_help = (snap.mach[iL]**2+1)/(snap.mach[iL]**2-1)

            Hmach   = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask, HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.mach[iL])
            Halpha  = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask, HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*alpha_help)
            Hrho_up = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask, HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.rup[iL]*fac_rho)
            Htemp   = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask, HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.u[iL]*fac_T)
            Harea   = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask, HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.area[iL])
            Hmag    = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask, HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radi[iL]*x.B[iL])
            Hpre    = SPH_binning(snap, posrot, dinf, iL, immask=DoG_mask, HSize_z=hsize, nbins=nbins, weights=lambda x: s_radio_SI*x.radiPre[iL])

            allflux = np.asarray([])
            for relic in relics:
                relic.wMach  = Hmach[relic.pmask]
                relic.wT     = Htemp[relic.pmask]
                relic.wArea  = Harea[relic.pmask]
                relic.wAlpha = Halpha[relic.pmask]
                relic.wB     = Hmag[relic.pmask]
                relic.wPre   = Hpre[relic.pmask]
                relic.wRho_up = Hrho_up[relic.pmask]
#                    relic.wRho        =  Hrho [relic.pmask]
#                    relic.wRho_down   =  Hrho_down[relic.pmask]  
#                    relic.wT_up       =  Htemp_up[relic.pmask]  
#                    relic.wT_down     =  Htemp_down[relic.pmask]  

                relic.wDoG_rel = DoG_rel[relic.pmask]
                allflux = np.concatenate((relic.sparseW, allflux), axis=0)
                relic.averages_quantities()

            """Save maps"""
            allflux = allflux.flatten()

            """ I couldn't come up with something better to take the inverse """
            mask = np.ones(snap.rho.shape, dtype=bool)
            mask[iL] = 0
            Subtracted, xedges, yedges = np.histogram2d(-posrot[mask,0], -posrot[mask,1], weights=s_radio_SI*snap.radi[mask],
                                                        range=[[-hsize,hsize], [-hsize,hsize]], bins=nbins)
            Subtracted += H1*(1-DoG_mask)
            Subtracted_conv = dinf.convolve_map(Subtracted)
            if extract_subtracted:
                relics_subtracted = relex.RelicExtraction(Subtracted_conv, z, GCl=gcl, dinfo=dinf,
                                                          rinfo=cbclass.RelicRegion('', [], rtype=1))  # , faintexcl=0.4, Mach=Hmach, Dens=Hdens,
                for relic in relics_subtracted:
                    relic.wMach = Hmach[relic.pmask]
                    relic.wT = Htemp[relic.pmask]
                    relic.wArea = Harea[relic.pmask]
                    relic.wAlpha = Halpha[relic.pmask]
                    relic.wB = Hmag[relic.pmask]
                    relic.wPre = Hpre[relic.pmask]
                    relic.wRho_up = Hrho_up[relic.pmask]
                    #                    relic.wRho        =  Hrho [relic.pmask]
                    #                    relic.wRho_down   =  Hrho_down[relic.pmask]
                    #                    relic.wT_up       =  Htemp_up[relic.pmask]
                    #                    relic.wT_down     =  Htemp_down[relic.pmask]

                    relic.wDoG_rel = DoG_rel[relic.pmask]
                    relic.averages_quantities()
                gcl.compacts = relics_subtracted

            if saveFITS:
                """ Here the maps are already masked with the detection region """

                smt(task='WriteFits_[writes,sub]')
                if log: print('###==== Step 4:  Preparing FITS file & folders ====###')

                parlist = (savefolder, gcl.mockobs.z_snap * 1000, gcl.name, Rmodel.id)
                gcl.maps_update(H1, 'Raw', '%s/maps/z%04.f/%s-%04i_native.fits' % parlist)
                gcl.maps_update(IM1, 'Diffuse', '%s/maps/z%04.f/%s-%04i.fits' % parlist)
                gcl.maps_update(Subtracted, 'CompModell', '%s/maps/z%04.f/%s-%04i_compact.fits' % parlist)
                gcl.maps_update(Subtracted_conv, 'Subtracted',
                                '%s/maps/z%04.f/%s-%04i_compactObserved.fits' % parlist)
                if len(relics) > 0:
                    gcl.maps_update(gcl.Mask_Map(Hmach, normalize=allflux), 'Mach',
                                    '%s/maps/z%04.f/%s-%04i_mach.fits' % parlist)
                    gcl.maps_update(gcl.Mask_Map(Hrho_up, normalize=allflux), 'RhoUp',
                                    '%s/maps/z%04.f/%s-%04i_rhoup.fits' % parlist)
                    gcl.maps_update(gcl.Mask_Map(Htemp, normalize=allflux), 'Temp',
                                    '%s/maps/z%04.f/%s-%04i_temp.fits' % parlist)
                    gcl.maps_update(gcl.Mask_Map(Hmag, normalize=allflux), 'B',
                                    '%s/maps/z%04.f/%s-%04i_B.fits' % parlist)
                    gcl.maps_update(gcl.Mask_Map(Hpre, normalize=allflux), 'PreRatio',
                                    '%s/maps/z%04.f/%s-%04i_prerat.fits' % parlist)
                gcl.maps_update(DoG_rel, 'DoG_rel', '%s/maps/z%04.f/%s-%04i_DoG_rel.fits' % parlist)
                gcl.maps_update(DoG_mask, 'DoG_mask', '%s/maps/z%04.f/%s-%04i_DoG_mask.fits' % parlist)

            """ PhD feature --> plot the DoF images in a subplot
##                import matplotlib.pyplot as plt   
##                with np.errstate(divide='ignore', invalid='ignore'):
##                    DoG_rel            = np.divide(np.abs(convolved_sigma2-convolved_sigma1)+1e-20,convolved_sigma2+1e-20)
##                pixR200 =        gcl.R200()/(gcl.cosmoPS*dinf.spixel)     
##                bou     =  gcl.R200()*1.5   # pixR200*1
##                f, ((ax1, ax2, ax3)) = plt.subplots(1, 3, figsize=(30,12)) #, sharey='row', sharex='col', sharey='row'
##                ax1.imshow( np.power((np.abs(convolved_sigma1-convolved_sigma2)), 1/1 )     , extent=(-hsize, hsize, -hsize, hsize)) #-bou+cen
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
#                """     """
            
            gcl.add_relics(relics) 
            PhD feature end """
    if writeClusters:
        """This is here because some outputs get lost in a multiprocessing heavy input/output queue process"""
        for gcl in GClrealisations_used:
            filename = 'GCl-%05i' % (gcl.mockobs.id)
            iom.pickleObject( (gcl, Rmodel), savefolder + '/pickled/', filename, append=False)


    if log: print( 'Finished with all efficency values')
    return True, smt, GClrealisations_used, Rmodel
