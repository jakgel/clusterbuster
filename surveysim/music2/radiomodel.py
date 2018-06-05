#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 10:52:22 2017

@author: jakobg
"""

from __future__ import division,print_function

import surveysim.music2.interpolate as interpolate
import surveysim.music2.loadsnap    as loadsnap                 

import numpy                     as np
import clusterbuster.mathut      as math
import clusterbuster.iout.misc   as io
import clusterbuster.surveyclasses  as cbclass

from scipy.special import expit

#from   bisect import bisect_left


#=== main functions
def PrepareRadioCube(snap, psiFile='Hoeft_radio/mach_psi_table.txt', machFile='Hoeft_radio/q_mach_machr_table.txt', log=False):
    ''' 
    This does some preparation steps that should be same, no matter, what model is used.
    I strongly recommend to use interpolated files like Hoeft_radio/mach_psi_tablefine(10,3).txt & Hoeft_radio/q_mach_machr_tablefine(10,3).txt
    - finner interpolation steps will slow down the computation.'''
    
    smt = io.SmartTiming(rate=5e4); smt(task='PrepareRadioCube')

    if log: print('###==== Step 0a:  Prepare cluster data cubes ====###')
    #==== Load psiFile   for psi factor & machfile for mach-numbers conversion factors
    H_mach      = np.loadtxt(machFile,skiprows=0) 
    H_psi       = np.loadtxt(psiFile,skiprows=0)[:,1::] # you wont get the temperature values ... this is why we read them separetely

    psi_x, psi_y = interpolate.LoadFile_psi(psiFile)
    
    # First: Apply a filter for points with mach number > psi_y[0] = 1.23 equivalent to the minimal value in the machlist
    # I implement a lower cut of m=1.5 just because it decreases the number of computed particles quite significantly
    MF   = np.where(snap.mach > 1.5)
    
    #==== Gets conversion factors
    rho_to_ne = loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb()  # [Msol parsec-3] com+h--> [electrons cm-3] physical
    U_to_keV  = loadsnap.conversion_fact_gadget_U_to_keV(  snap.head )   
    
    # Derive Temperatur (in subarray)
    T         = U_to_keV*snap.udow[MF]              # in [keV]  
 
    #==== Finding the closest corresponding value to table entries --> x is for Temperature, y stands for the mach number
    results_x = math.find_closest(psi_x, T                     )
    results_y = math.find_closest(psi_y, snap.mach[MF]*1.045   ) #

    s         = H_mach[results_y,4]         # Electron energy distribution spectral index as function of mach number    
     
     
    # ... get an idea of the corresponding iluminated area 
    f         = 6.5                                               # Shock area fraction due to tab. 1 in  Hoeft+2008,  also consider M -> M *1.045 due to  Hoeft+2008
    h         = snap.hsml[MF]*1e-3                                # com+h [kpc] --> com+h [Mpc]  Hydrodynamical smoothing length, i.e. Size of kernel, determined as a measure of density  IMPORTANT: The term h^-1*expansion factor (because of comoving frame) is added later in CreateRadioCube!
    N_k       = 64.                                               # The exact number of particles in the smoothing kernel is yet not known by me --> ask e.g. Sebastian Nuza
    factor_A  = loadsnap.comH_to_phys( snap.head )                      # Smoothing kernel part B; 0.7 steems from h=0.7; 1/(1+z) is the expansion factor; second part of computing A_i   
    A_i       = f*h*h/N_k*factor_A**2                         # in [Mpc]^2 ... Inspired by equ 18 Hoeft+2008 (http://mnras.oxfordjournaloadsnap.org/content/391/4/1511.full.pdf)                            
    
    
#    ''' DEBUGGING to infer the internal smoothing length and to compare it with other equations '''
#    for luck in zip(snap.hsml[::30000], snap.rho[::30000]):
#             print(luck[0]*factor_A,luck[1]*rho_to_ne*1.e4)
             
    # Get other values for radio emission formula
    rho_e     =  snap.rdow[MF]*rho_to_ne*1.e4 # in [electrons 10^-4 cm^-3] az z=0
    Xi_e      =  1.                           # Xi_e : energy fraction of suprathermal electrons, around 10^(-6...-4)  by default set to one and rescalled later on
    
    ''' If the Gelszinnis model if PREs is not used, this can be omitted '''
    snap.DSAPsi = snap.rdow*0 
    snap.DSAPsi = H_psi[results_y,results_x]
    '''==='''
    
    #=== Compute Radio Emission: Using Eq. 32 from Hoeft-Brueggen 2007
    snap.radi = np.array(snap.rdow*0,dtype='float64')    #Used to initiate a numpy array of the same size and shape as the others; the large float currently comes from radio emision which is to high
    if log: print("PrepareRadioCube (np.min(A_i),np.max(rho_e),np.max((Xi_e/0.05)),np.max((T/7)),np.max(H_psi[results_y,results_x]): %10.2e %10.2e %10.2e %10.2e %10.2e" % (np.min(A_i),np.max(rho_e),np.max((Xi_e/0.05)),np.max((T/7)),np.max(H_psi[results_y,results_x])))
#    print(':::', np.sum(A_i), np.sum(rho_e), np.sum(Xi_e), np.sum(np.power(T/7.,1.5)), np.max(H_psi[results_y,results_x]) )
    snap.radi[MF] =  6.4e34 * A_i * rho_e * (Xi_e/0.05) * np.power(T/7.,1.5)*snap.DSAPsi     # in  [erg/s/Hz] divided by (factur_nu*factor_B)
#    print(':_:_:', np.sum(snap.radi[MF]) )
    
    # The spectral index of the electron energies s / Equals two times the injection index
    snap.s     = snap.rdow*0 
    snap.s[MF] = s
   
    #=== Add the area
    snap.area      = np.array(snap.rdow*0,dtype='float64')   #Used to initiate a numpy array of the same size and shape as the others; the large float currently comes from radio emision which is to high
    snap.area[MF]  = A_i  
          
    
    # The spectral index of the (down stream integrated) radio emission is (1 - s)/2, SO THAT ALPHA FOR RADIO RELICS SHOULD BE NEGATIVE
    #snap.alpha =   snap.rdow*0
    #snap.alpha  = -(s-1.)/2.

    return (snap, MF)  #better a suparray like snap.radi and snap.alpha
    
    

def PiggyBagSnap(snap, extended=True):
    ''' Chooses just a few properties of the snap to reduce the load of interprocess communication  later on 
    
    extented = True: yields additional information for deeper analysis
    
    '''

    cutsnap          = loadsnap.Snapshot(snap.name) 
    cutsnap.head     = snap.head
    cutsnap.rup      = snap.rup
    cutsnap.rdow     = snap.rdow
    cutsnap.mach     = snap.mach 
    cutsnap.s        = snap.s
    cutsnap.radi     = snap.radi
    cutsnap.pos      = snap.pos
      
    if extended:
        cutsnap.udow   = snap.udow
#        cutsnap.rdow   = snap.rdow
        cutsnap.mach   = snap.mach
        cutsnap.area   = snap.area
        
        '''DEVELOPMENT I --> this is needed, because it PiggyBagSnap is also used for Cutsnap in run survey, which leads to PiggyBagSnap_cut'''
        cutsnap.hsml   = snap.hsml
        cutsnap.uup    = snap.uup
        
        '''DEVELOPMENT II  Because in CreateMockOps, cuts are set due to the density and temperature of the projectile'''
        cutsnap.rho    = snap.rho
        cutsnap.u      = snap.u
        
        '''DEVELOPMENT III
        Because we create a model with PREs, where the Gelszinnis model uses a modification of this factor. If not used, this can be ommited'''
        cutsnap.DSAPsi = snap.DSAPsi    
 
    return cutsnap

def PiggyBagSnap_cut( snap, procsnap, cutradio=None, cutmask=None): 
    ''' Cuts snap to relevant particles and parameters based on a flat cut in snap 
        Since 2018, we apply a flat cut based on the Mach-number of the particles.
    '''

    if snap.rdow.shape != procsnap.radi.shape:
        print('Something went wrong with the array! ..., snap.shape, procsnap.shape', snap.rdow.shape, procsnap.radi.shape)

    # Get the mask
    if cutmask is None:
        cutmask = np.where( procsnap.radi > cutradio)
        
    # Pick the relvant data
    cutsnap            = loadsnap.Snapshot(snap.name) 
    cutsnap.name       = snap.name
    cutsnap.head       = snap.head
    cutsnap.rup        = snap.rup[cutmask]
    cutsnap.rdow       = snap.rdow[cutmask]
    cutsnap.mach       = snap.mach[cutmask] 
    cutsnap.s          = snap.s[cutmask]
    cutsnap.radi       = snap.radi[cutmask]
    cutsnap.hsml       = snap.hsml[cutmask]
    cutsnap.uup        = snap.uup[cutmask]
    cutsnap.udow       = snap.udow[cutmask]
    cutsnap.pos        = snap.pos[cutmask]
    
    '''DEVELOPMENT II  Because in CreateMockOps, cuts are set due to the density and temperature of the projectile'''
    cutsnap.rho    = snap.rho[cutmask]
    cutsnap.u      = snap.u[cutmask]

    return cutsnap
  
  

def CreateRadioCube( snapred, Rmodel, z,  nuobs=1.4, logging = True):
    ''' The returned radio luminousity is the radio luminousity in the restframe; Compare Nuza+ 2012 Equ (1)
        Also see Section 5 Araya-Melo 2012 et al.  
     Please mind that the density (-->B) estimate is set to the redshift of the cluster snap, but Bcmb and observing frequency to the redshift of the mockobs
     
     
     
    '''
    compress = Rmodel.compress
    kappa    = Rmodel.kappa
    B0       = Rmodel.B0
    
    snap, MF = snapred
    smt = io.SmartTiming(rate=5e4); smt(task='CreateRadioCube_[sub]')
    
    
    # Get the density right
    nurest     = float(nuobs)*(1+z)  # [GHz] It differs from observing frequency!
    z_snap     = 1/snap.head['aexpan'] -1
    f_cor      = (z+1.)/(z_snap+1.)   # Now we compute the cube at the correct redshift ... we just have to compute a correction factor! :-)

    rho_conv    = f_cor**3*loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb()*1e4 # in [electrons 10^-4 cm^-3] az z=z_obs
    T_conv      = loadsnap.conversion_fact_gadget_U_to_keV(  snap.head )/8.61732814974056e-08     # K
    # Compute magnetic field
    
    if compress == -1: #Uses the scaling by Nuza+2016
        R          = 1
        rho_e      = rho_conv*snap.rdow[MF] #Downstream density
        B          = B0*np.power(rho_e, kappa)  # in [muG] - using formula for B-field strength
    else:
        R          = np.divide(snap.rdow[MF],snap.rup[MF])**compress  # Compress of 1 would show the same behavior like the Nuza model, just for an lower magnetic field
        rho_e      = rho_conv*snap.rup[MF]        # Upstream density

        ''' DEVELOPMENT '''
        if Rmodel.B_para == 'press': 
            B  = B0*np.power(rho_e*snap.u[MF]/1e6, kappa)*R  # in [muG] - using formula for B-field strength  and compression
        elif Rmodel.B_para == 'dens': 
            ''' Nuza+2017 parametrisation '''
            B  = B0*np.power(rho_e, kappa)*R  # in [muG] - using formula for B-field strength  and compression
        else:
            print('Parametrization of the magnetic field is unknown. Choose BFieldModel.B_para.') 
        B[B > 1e5] = 1e5                          # set a maximumB field, this has to be done as in some cases for odd magnetic values, we get overflow of magnetic fields valus, as it seems .... 
    
    snap.B      = np.array(snap.u*0,dtype='float64')   #Used to initiate a numpy array of the same size and shape as the others; the large float currently comes from radio emision which is to high
    snap.B[MF]  = B
          
    Bcmb      = 3.25*(1+z)**2.                    # in [muG] - taken from equ. shortly before 2.43 in  Galactic and Intergalactic Magnetic Fields by Ulrich Klein, Andrew Fletcher - 2014 , Science]
    # Compute additional factors needed for radio flux computations
    factor_nu  = np.power( (nurest/1.4) , (-snap.s[MF]/2.) ) # GHz
    factor_B   = np.power(B, 1.+snap.s[MF]/2.) / ( B**2. + Bcmb**2. )       # Magnetic field strength 
    factor_fudge =   (f_cor)**3/(f_cor)**2 / (1+z)**0.0    #3.5             # First time for density factor, second term for area third term is just to make it fit ... It is set to 1 (no correction) as I interpret the radio emissivity in sebastians cube as the one for the observers frequency

    f_boost     = np.array(snap.rdow*0,dtype='float64')     
    '''Is based on the idea of a boosting-factor, and is related to more recent discussion of our working group.
       The presumption is that an additional CR population exists. In the following model it is assumed that
       it steems from accretion shocks. M.Hoeft deriffed a boosting factor 'f_boosed' based on this assumption.
       
       remark:
       Some authours say, that merger shocks put more energy (and even more highly relativistic particles) in the cluster. This model is not considered here
    '''
    if Rmodel.pre:
        
        if isinstance(Rmodel, cbclass.PreModel_Hoeft):
        
            delta                  = np.log10(rho_e/Rmodel.n0)/np.log10(Rmodel.n1/Rmodel.n0)
            taccr                  = delta*Rmodel.t1+(delta-1)*Rmodel.t0
            taccr[rho_e<Rmodel.n0] = Rmodel.t0
            taccr[rho_e>Rmodel.n1] = Rmodel.t1
            gamma_boost            = 2.4e4/taccr/( (B/R)**2. + Bcmb**2.)  #Magnetic field before enhancement through compression
            f_boost[MF]            = Rmodel.ratio*gamma_boost**(snap.s[MF]-2)
        
        elif isinstance(Rmodel, cbclass.PreModel_Gelszinnis):

            '''modify xhi,
                 #add a certain fraction (due to amount of preexisting electrons) to emission
                 #Add exponential term to emissivity of xhi --> strengthens the emission of low mach number shocks
                 #Also (log?)normal normalization of plasma --> should saturate at a certain fraction
                 #Influence for high mach number shocks should roughly add the same amount of emission as the thermal pool in average.
            '''
            f_expid = expit(  (Rmodel.sigmoid_0-np.log10(snap.rdow[MF]*snap.u[MF]*(rho_conv*T_conv)))/Rmodel.sigmoid_width  )
            s       = 1     
            if Rmodel.p_sigma  > 0:
                s   = np.random.normal(0, Rmodel.p_sigma)                                                
            f_boost[MF] = s*np.power(snap.DSAPsi,Rmodel.PREexpo-1)*f_expid
        else:
            print('The model of pre-existing electrons is unknown. Choose the inheriting class of ObjectClasses::Rmodel.')

                     
    #=== Compute Radio Emission: Using Eq. 32 from Hoeft-Brueggen 2007, part B
    ergSI         =  1.e-7                                                     # [erg/s per W]
    snap.radi[MF] =  ergSI * snap.radi[MF]*factor_nu*factor_B*factor_fudge     # in  [W/Hz]
    snap.radiPre  =  snap.radi*f_boost                                         # in  [W/Hz]
    #print 'snap.s[MF]', snap.s[MF] DEBUGGING CHEAT here!
    if logging: print('  Total radio power in cube: %5.2e W/Hz at %5.2f GHz rest frame frequency' %  ( np.sum(snap.radi[MF]), nurest))
    radiofolder = 'ToImplement'  #strSn.replace('shocks','radio')  
    
    return (snap, smt, radiofolder)  # A tuple, only [0:2] is relevant
    


def PlotParticleFilter(snap, iL, savefile, particleW=lambda x: x.mach**2, labels=False,summed=True):
    ''' Plots filtered and unfiltered particle due to a given mask iL of a snapshot and 
        compares it to another weight
        
        It was added to facilitate a reasonable cut of particles in the phase space (density- u)
    '''
    fac   = loadsnap.conversion_fact_gadget_rho_to_nb( snap.head )*loadsnap.conversion_fact_ne_per_nb()
    facU  = loadsnap.conversion_fact_gadget_U_to_keV(  snap.head )/8.61732814974056e-08  # to K
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    

    ''' I couldnt come up with omething better to take the inverse '''
    mask = np.ones(snap.rho.shape,dtype=bool)
    mask[iL]=0
    
    fig = plt.figure(figsize=(7, 3))
    ax1 = fig.add_subplot(121, title='particles')
    
    plotrange = [[-10,2], [2.0,10]]
    bins      = (90,60)  #(30,30)
    
    '''Image 1 excluded'''
    H_mass_full, xedges, yedges = np.histogram2d ( np.log10(snap.rho*fac), np.log10(snap.u*facU), range=plotrange, bins=bins) #weights=np.ones( (iL.shape[0])), 
    H_mass_full = H_mass_full.T  # Let each row list bins with common y range.
    plt.imshow(H_mass_full, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], norm=mpl.colors.LogNorm(),alpha=0.65)
    ax1.set_ylabel('$\\log_{10}\\left(T/\mathrm{K}\\right)$')
    ax1.set_xlabel('$\\log_{10}\\left(n_e/\mathrm{cm^{3}}\\right)$')
    
    '''Image 1 '''
    H_mass, xedges, yedges = np.histogram2d ( np.log10(snap.rho[iL]*fac), np.log10(snap.u[iL]*facU), range=plotrange, bins=bins) #weights=np.ones( (iL.shape[0])), 
    H_mass = H_mass.T  # Let each row list bins with common y range.
    plt.imshow(H_mass, interpolation='nearest', origin='low'    ,extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], norm=mpl.colors.LogNorm())
    

    ax2 = fig.add_subplot(122, title='particles weighted')
    '''Image 2'''
    H_weighted, xedges, yedges = np.histogram2d ( np.log10(snap.rho*fac), np.log10(snap.u*facU), weights=particleW(snap), range=plotrange, bins=bins) #weights=np.ones( (iL.shape[0])), 
    H_weighted = H_weighted.T  # Let each row list bins with common y range.

    kwargs = {}
    if labels: 
        kwargs = {'norm':mpl.colors.LogNorm()}

    plt.imshow(H_weighted, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], **kwargs) #, norm=mpl.colors.LogNorm()
#    
#    H_mass, xedges, yedges = np.histogram2d ( np.log10(snap.rho[iL]*fac), np.log10(snap.u[iL]*facU), weights=particleW(snap)[iL], range=plotrange, bins=(30,30)) #weights=np.ones( (iL.shape[0])), 
#    H_mass = H_mass.T  # Let each row list bins with common y range.
#    plt.imshow(H_mass, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]) #, norm=mpl.colors.LogNorm()
#    
#    H_mass, xedges, yedges = np.histogram2d ( np.log10(snap.rho[mask]*fac), np.log10(snap.u[mask]*facU), weights=particleW(snap)[mask], range=plotrange, bins=(30,30)) #weights=np.ones( (iL.shape[0])), 
#    H_mass = H_mass.T  # Let each row list bins with common y range.    
#    plt.imshow(H_mass, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],alpha=0.65) #, norm=mpl.colors.LogNorm()
   
    ax2.set_xlabel('$\\log_{10}\\left(n_e/\mathrm{cm^{3}}\\right)$')

    if labels:
        ax2.set_xlabel('$\\log_{10}\\left(\\rho/\mathrm{cm^{-3}}\\right)$')
#        
#        import matplotlib
#        from matplotlib.patches import Polygon
#        from matplotlib.collections import PatchCollection
#        
#        fig, ax = plt.subplots()
#        patches = []
#        num_polygons = 2
#        num_sides = 5
#        
#        for i in range(num_polygons):
#            polygon = Polygon(np.random.rand(num_sides ,2), True)
#            patches.append(polygon)
#        
#        p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
#        
#        colors = 100*np.random.rand(len(patches))
#        p.set_array(np.array(colors))
#        
#        ax2.add_collection(p)
        ax2.text(0.22, 0.32,'virial', fontweight='bold', horizontalalignment='center',
              verticalalignment='center', color='white',rotation=30,
              transform=ax2.transAxes)
   
        ax2.text(0.43, 0.65,'hot ICM', fontweight='bold', horizontalalignment='center',
              verticalalignment='center', color='white',rotation=30,
              transform=ax2.transAxes)
    
        ax2.text(0.76, 0.32,'cooled', fontweight='bold', horizontalalignment='center',
              verticalalignment='center', color='white',rotation=30,
              transform=ax2.transAxes)
    
    plt.savefig(savefile)
    
    
    outfile = '/data/ClusterBuster-Output/Test'
    
    
    
    if summed:
               
        H_mach, xedges, yedges = np.histogram2d ( np.log10(snap.rho*fac), np.log10(snap.u*facU), weights=snap.mach, range=plotrange, bins=bins) #weights=np.ones( (iL.shape[0])), 
        H_mach = H_mach.T  # Let each row list bins with common y range.
           
        print(outfile)   

        import time
        attempts = 0
        while attempts < 10:
            try:
                Hadd          = np.load(outfile+'_particles.npy')
                Hadd_mach     = np.load(outfile+'_mach.npy')
                Hadd_weighted = np.load(outfile+'_weighted.npy')
                print('__________',np.sum(Hadd),np.sum(Hadd_weighted))
                break
            except:
                attempts += 1
                Hadd          = np.zeros_like(H_mass)
                Hadd_mach     = np.zeros_like(H_mass)
                Hadd_weighted = np.zeros_like(H_mass) 
                time.sleep(0.2)

                             
        Hadd         = H_mass_full +Hadd   
        Hadd_mach    = H_mach      +Hadd_mach   
        Hadd_weighted= H_weighted  +Hadd_weighted   
        np.save(outfile+'_particles.npy', Hadd) 
        np.save(outfile+'_mach.npy' , Hadd_mach)   
        np.save(outfile+'_weighted.npy' , Hadd_weighted)         
  
        
        x  = np.linspace(-10, 11, num=500, endpoint=True)
#        x1 = np.linspace( -5, -3, num=500, endpoint=True)    
#        x2 = np.linspace(-10, -5, num=500, endpoint=True)
#        x3 = np.linspace( -3, 11, num=500, endpoint=True)
        
        y1 =   3.3 - 0.65*x
        y2 = -11.0 - 3.5*x
        y3 =   6.6 + 0.5*x
        y  = np.maximum(np.maximum(y1,y2),y3)
        
        '''Image 1 excluded'''
        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(131, title='a)') #particles
        plt.imshow(Hadd, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], norm=mpl.colors.LogNorm())
        ax1.set_ylabel('$\\log_{10}\\left(T/\mathrm{K}\\right)$')
        ax1.set_xlabel('$\\log_{10}\\left(n_e/\mathrm{cm^{3}}\\right)$')


        ax2 = fig.add_subplot(132) #, title='average Mach'
        im2 = plt.imshow(Hadd_mach/Hadd, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], vmax= 50, norm=mpl.colors.LogNorm())
        ax2.set_xlabel('$\\log_{10}\\left(n_e/\mathrm{cm^{3}}\\right)$')
#        divider = make_axes_locatable(ax2)
#        cax2 = divider.new_vertical(size="5%", pad=0.7, pack_start=True)
#        cax2 = fig.add_axes([0.42, 0.8, 0.2, 0.03]) 
#        cb2  = fig.colorbar(im2, format='%.0f', ticks=[3, 10, 30], cax = cax2, orientation="horizontal")  #label='average Mach', 
#        frame = plt.gca()
        ax2.yaxis.set_ticklabels([])

        ax3 = fig.add_subplot(133, title='c)') #summed flux
        kwargs = {}
        if labels: 
            kwargs = {'norm':mpl.colors.LogNorm()}
#        print('__',np.mean(np.log10(Hadd_weighted)))
        plt.imshow(Hadd_weighted, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], norm=mpl.colors.LogNorm(), vmin=1e25, vmax=2e30, **kwargs) #, norm=mpl.colors.LogNorm()
        ax3.set_xlabel('$\\log_{10}\\left(n_e/\mathrm{cm^{-3}}\\right)$')
        ax3.yaxis.tick_right()
        
        
        ax1.autoscale(False)
        ax2.autoscale(False)
        ax3.autoscale(False)  
#        ax1.plot(x1,y1,c='r')
#        ax1.plot(x2,y2,c='r')
#        ax1.plot(x3,y3,c='r') 
        ax1.plot(x,y,c='r') 
        ax2.plot(x,y,c='r') 
        ax3.plot(x,y,c='r') 
        
        ax2.text(0.22, 0.32,'virial', fontweight='bold', horizontalalignment='center',
              verticalalignment='center', color='white',rotation=30,
              transform=ax2.transAxes)
   
        ax2.text(0.43, 0.65,'hot ICM', fontweight='bold', horizontalalignment='center',
              verticalalignment='center', color='white',rotation=30,
              transform=ax2.transAxes)
    
        ax2.text(0.76, 0.32,'cooled', fontweight='bold', horizontalalignment='center',
              verticalalignment='center', color='white',rotation=30,
              transform=ax2.transAxes)
        
        plt.savefig('/data/ClusterBuster-Output/Test_SummedSurvey.pdf')
        
        
        
        
'''DEVELOPMENT PHD'''    

if __name__ == "__main__":
    import clusterbuster.IOutil                    as iout
    snap    = iout.unpickleObject('/data2/MUSIC-2/snaps/SHOCKS_00001/clusterSUBSET.00001.snap.014.shocks')   
    iL      = np.where(snap.rup > 0)
    PlotParticleFilter(snap, iL, '', summed=True)            
                
'''          '''