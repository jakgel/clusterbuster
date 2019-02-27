#!/usr/bin/env python

"""
Created on Mon Dec 11 17:34:32 2017
@author: jakobg
This module includes custom classes used to create an logical-hierarchy of the properties of galaxy clusters, the objects within and 
the means by which they were observed: Surveys, Observation configurations, and synthesic observations settings.

These object classes are pickled. An fully completed project would include that the survey object would be throuroghly parset into a database.
"""

from __future__ import division,print_function,absolute_import

import operator
import os
import copy
import pandas as pd  
#import warnings

import ephem as ephem
import numpy as np
import clusterbuster.mathut as maut
import clusterbuster.surveyut as surut
import clusterbuster.fitsut as fitsut
import clusterbuster.dbclasses as dbc

import scipy.stats.mstats as mstats #.gmean
import scipy.ndimage as ndi
import matplotlib.patches as patches
import NFW.mass_concentration as massC

#from .mathut  import *

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM   # from astropy.cosmology import WMAP9 as cosmo
from scipy import sparse
from scipy.ndimage.filters import gaussian_filter1d


def get_truth(inp, relate, cut):
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '=': operator.eq}
    return ops[relate](inp, cut)


def replaceNone(variable, default):
    
    if variable is None:
        return default
    else:
        return variable

class Survey(object):
    """
    Simple class for representing a  survey for radio relics in galaxyclusrers galaxy clusters
    It provides functions and methodologies to run ...
    To initialize it needs a  list of galaxy clusters to be assigned.
    ... might add astropy functionalities --> including using astropy quantities (quantity + unt)
    ... might write it for python 3.X
    """
    
    def __init__(self, GCls, survey, cnt_levels=[0], synonyms=[], Rmodel=None, emi_max=False, scatterkwargs={}, histokwargs={},
                 saveFITS=True, savewodetect=False, dinfo=None, hist_main=None, surshort='surshort', logfolder='', outfolder = None):
        """
        parameters
        """
        self.name      = survey
        self.name_short= surshort
        self.GCls      = GCls      # A list of galaxy clusters
        self.DetInfo   = DetInfo   # A specified DetInfo as proxy for an homogeneous survey; Detinfo include a beam, and a noiseMAP --> properties should be taken from the neirest entry
        self.dinfo     = dinfo
        self.hist_main = hist_main  # Histogram of all binned objects/ without any initialization, works as an proxy

        self.filteredClusters = None
        
        """This should be handled in the case of a mock survey"""
        """So there should be a nock survey and a real survey, both as childs from a generall survey"""
        self.Rmodel = replaceNone(Rmodel, RModel([], simu=False))   # A B-field model
        
        # Total measures, you can also write an function for that

        self.P_pro  = 0        # Pro ratio
        self.P_anti = 0        # Pro ratio
        self.F_pro  = 0        # Pro ratio
        self.F_anti = 0        # Pro ratio

        # Output related properties
        if outfolder is None: outfolder = '/data/ClusterBuster-Output/%s' % (survey)
        self.outfolder = outfolder
        self.logfolder = logfolder
        self.saveFITS = saveFITS
        self.synonyms = synonyms
        self.savewodetect = savewodetect
        
        # plot Emission I: Contours & colorscale
        self.cnt_levels = cnt_levels
        self.emi_max    = emi_max      # Maximal scale of emission
        self.m_cnt      = 2            # Multiplier for contour levels
        
        # plot Relics & Clusters: Scatterplots & Histogram
        self.scatter_backgr = False
        self.scatterkwargs  = scatterkwargs  # A list of kwargs for the scatter plots
        self.histokwargs    = histokwargs    # A list of kwargs for the histogram plots
        
        # Plot olar
        self.expScale   = 0.75
        self.expA       = 1.
        self.seed_dropout = None
        
    def set_binning(self, histo=None):

        if histo is None:
            histo = self.histo
        for GCl in self.GCls:
            GCl.histo = histo
        
    def polar(self, aligned=True, normalize=True, minrel=1,  positive=False,  mirrored=False, mode='flux', zborder=0,
              ztype='>', conv=0.127,  **kwargs):

        """ 
        Input: needs to histogram of the survey to be set, radial axis: 2 N, polar axis: 4 M
        conv: sigma of convolution of the radial bin in R200
        
        Returns the weighted polar Histogramm of the survey, currently in Output class 
        
        set normalize=True: if you want to get the cumulative signal divided by the number of relic hosting clusters
        
        Output an Tuple of 4 outputs:
          halfHist_new (np.array) : halfHist for analysis issues, shape stays same
          radial & ticks_r          : (radial,[-h for h in reversed(Histo.ticks[1])]+[h for h in Histo.ticks[1]]),
          halfHist_plot           : halfHist for plotting purposses, shape may change depending on the output format
          sigstats (float)       :  total accumulated signal
          mesh                    :  total accumulated signal along

        """
          
        """DEVELOPMENT"""
        ndist = 35  #-->deprecated! but a better umber (nicer looking then the raw grid)
        distmax = 3500  # very important, should somehow be part of the histogramm class
        # This is a design desicion ... in the standart case, the histogramm should ourient itself on the histogram variable of the Survey class
        mesh = np.zeros_like(self.hist_main.hist.T)
        """DEVELOPMENT END"""
        
        hist_main = np.zeros_like(self.hist_main.hist.T)
        sigstats = []
        
        if self.filteredClusters is None:
            self.FilterCluster(minrel=minrel, zborder=zborder, ztype=ztype)
        for ii, GCl in enumerate(self.filteredClusters):
            GCl.histo = self.hist_main
            GCl.updateInformation(Filter=True)
            if GCl.histo is not None and np.sum(GCl.histo.hist) > 0:
    
                Histo = GCl.histo
                hist_shift = np.roll(Histo.hist.T, -int(aligned*(GCl.relic_pro_index)), axis=1)  # This was a bug:/ AreaHist**(self.expA)
                """ The scale is very important, it would be good to normalize the flux we do this by implementing the expScale parameter.
                This is flexible enough to give different weights to bright/faint relics. Another question for statistical analysis is, if
                the average  or an integrated value is more interesting.
                """  
                
                if mode == 'flux':
                    scale = 1. 
                elif mode == 'power':
                    scale = 1. * GCl.flux2power_woRedening()    
                
                signal = np.divide(hist_shift*scale, (np.sum(hist_shift)*scale)**self.expScale)
                hist_main += signal
                bin_cluster = int(GCl.R200()/distmax *ndist)
                mesh[:, bin_cluster] += np.sum(signal, axis=1) 
                sigstats.append(np.sum(signal))
                
        if len(sigstats) > 0:
            if normalize: 
                hist_main = hist_main/max(len(sigstats), 0)
            halfHist = np.add(hist_main[:,0:int(hist_shift.shape[1]*1/2)], np.fliplr(hist_main[:,int(hist_shift.shape[1]*1/2):]) )

            """ Transformation for plotting: Only works if histo was already set"""
            if mirrored:
                halfHist_plot = np.add(np.fliplr(hist_main[:,0:int(hist_shift.shape[1]*1/2)]), hist_main[:,int(hist_shift.shape[1]*1/2):] )                  
            else:
                halfHist_plot = np.copy(halfHist)
                  
            # We compute a weighting array to weight to bins at larger radii, which are larger down
            inner = Histo.bins[1][0:-1]
            outer = Histo.bins[1][1::]
            angle = Histo.ticks[0]
            angles, z_inner = np.meshgrid(angle, inner, sparse=True)
            angles, z_outer = np.meshgrid(angle, outer, sparse=True)
            AreaHist = 2*np.pi*(2*np.pi)/len(angle)*(z_outer**2-z_inner**2)     
            halfHist_plot /= AreaHist**self.expA
            
            """ radially binned """
            part1 = np.sum(halfHist[:, int(halfHist.shape[1]/2):int(halfHist.shape[1]+1)], axis=1)[::-1]
            part2 = np.sum(halfHist[:, 0:int(halfHist.shape[1]/2+1)], axis=1)
            if positive:
                radial = np.sum( (part1, part2), axis=0)  # I would like to have np.flip
                ticks_r = [h for h in Histo.ticks[1]]
            else:
                radial = np.concatenate( (part1, part2), axis=0)  # I would like to have np.flip
                ticks_r = [-h for h in reversed(Histo.ticks[1])]+[h for h in Histo.ticks[1]]

            radial = gaussian_filter1d(radial, conv/self.hist_main.width[1], mode='constant') #1274
            #radial = radial/np.sum(radial)
            return halfHist, (radial, ticks_r), halfHist_plot, sigstats, mesh
        else:
            return None, (None, None), None, None, None
        
    def set_seed_dropout(self, seed_dropout=None):
        if seed_dropout is not None:
            self.seed_dropout = seed_dropout
        else:
            self.seed_dropout = np.random.RandomState()

    def FilterCluster(self, minrel=1, zborder=0, ztype='>', minimumLAS=0, GClflux=0, index=None, getindex=False,
                      verbose=False,  **kwargs):
        """ Gives all the cluster with the relics that fullfill given criteria """

        if index is not None:
            if index == 'All':
                return self.GCls
            else:
                return [self.GCls[i] for i in index]

        """ This part is implemented to give consistent results in plotting and metric etc. for a given seed_dropout """
        if self.seed_dropout is not None:
            state = self.seed_dropout.get_state()
            self.seed_dropout.set_state((state[0], state[1], 0, state[3], state[4]))

        GCls = [GCl.updateInformation(Filter=True) for GCl in self.GCls]

        results = [(ii,GCl) for ii, GCl in enumerate(GCls) if len(GCl.filterRelics(**kwargs)) >= minrel and
                   get_truth(GCl.z, ztype, zborder) and GCl.largestLAS() >= minimumLAS and GCl.flux() >= GClflux
                   and GCl.stoch_drop(self.seed_dropout)]

        """ Should also give an result, if no cluster fullfills the criterion """
        if len(results) > 0:
            indices, self.filteredClusters = map(list,zip(*results))
        else:
            indices, self.filteredClusters = [], []

        if verbose:
            print('____ # FilterCluster:', len(indices))

        if getindex:
            return indices
        else:
            return self.filteredClusters


    def fetch_totalRegions(self):
          
        return [GCl.relicRegions for GCl in self.GCls]
         
    def fetch_totalRelics(self, **kwargs):
        relics_list = []

        zborder = 0.05
        if self.filteredClusters is None:
            self.FilterCluster(zborder=zborder)

        for GCl in self.filteredClusters:
            relics_list += GCl.filterRelics(**kwargs)

        return relics_list
         
         
    def fetch_totalHisto(self):
          
        return sum([GCl.hist for GCl in self.GCls], 0)
     
    def fetch_pandas(survey, plotmeasures, logs=None,  surname=True, keys="label", vkwargs_FilterCluster={}, kwargs_FilterObjects={}):
        """ Return a panda array generated from the survey catalogue 
        
        survey: A ClusterBusterSurvey
        
        
        Examples for plotmeasures: 
            
        plotmeasures = [             
                        lambda x: x.alpha,
                        lambda x: dbc.measurand( x.Dproj_pix()/x.GCl.R200(), 'Dproj',label='$Dproj_rel$',  un = None )
                        ]
        """  

        if logs is None:
            logs = [True for measure in plotmeasures]  # a stub, you should be able to specify this
    
        list_full = []
        datalist = []
    
        """ This very ugly steps just creates a List of (relic,GCL) from the survey
        As in the future galaxy clusters will gain be a property of relics again this step will be shorter  and more readable"""
        for GCl in survey.FilterCluster(**vkwargs_FilterCluster):
                GCl.updateInformation()
                list_full.append([(GCl,relic) for relic in GCl.filterRelics()])
                
        list_full = [item for sublist in list_full for item in sublist]
            
        for GCl, relic in list_full:
    
            relic.asign_cluster(GCl)
            datavalues = []
            for measure, log in zip(plotmeasures, logs):
                data = measure(relic).value
                if log:
                    data = np.log10(data)
                datavalues.append(data)
            
            datalist.append(datavalues)

        """ Create a pandas dataframe """
        if keys == "label":
            columns = [measure(relic).labels(log=log) for measure,log in zip(plotmeasures, logs)]
        if keys == "dic":
            columns = [measure(relic).dic for measure, log in zip(plotmeasures, logs)]

        pdframe = pd.DataFrame(datalist, columns=columns)
        if surname:
            pdframe['Survey'] = survey.name

        if len(pdframe) < 3:
            print('Only %i elements in the pdframe, will skip this one' % (len(pdframe)))
            return 
        else:
            return pdframe


# Define a new object: relic
class Galaxycluster(object):           

    """
      Simple class for representing a galaxy cluster in our visible universe.
      ... Regions should contain relics, not other way around
      ... each relic (or region?) should a dinfo be assigned to
      ... region and dinfo might be outdated
      --- inconsistencies exist for the image arrays ... in the best case they would be converted to a new class to observation(s) (mockobs would then be a subset off obs)
      ... might add astropy functionalities --> including using astropy quantities (quantity + unt)
      ... might write it for python 3.X
      ... In add information: Roundish relics should e kept, with the notion that their are roundish
    """             

    def __init__(self, name, RA, Dec, z, M200=1e14, M500=0, M100=0, Mvir = 0, Lx=0, Lx_lit=0, flux_ps=0, flux_lit=0,
                 Prest100=0, relics=[], regions=[], halo=False, dinfo=None, ClassFlag=False, mockobs=None,
                                      Image=np.zeros((2,2)), status=None, Histo=None, reference=None, mapdic=dict()):
      
        """
        parameters
        """  
        
        self.name = name  #name.replace('_', ' ')
        self.RA   = dbc.measurand(RA, 'RA', un='deg')
        self.Dec  = dbc.measurand(Dec, 'Dec', un='deg')
        self.sky  = SkyCoord(self.RA(), self.Dec(), frame='icrs', unit='deg')
        self.z    = dbc.measurand(z, 'z',  label='$z$', un=None)  #Redshift

        """  We differentiate several (currently two) measures of the clsuter central position 
        One is the center of the dipole (we assume a double merger) of the cluster              sky_double
        The other one is the position of the X-Ray brightness peak of the more massive cluster  sky_main 
        """ 
        RA_d, Dec_d = RA, Dec # Curenty not implemented, DEBUGGING
        RA_m, Dec_m = RA, Dec # Curenty not implemented, DEBUGGING
        self.sky_double = SkyCoord(RA_d, Dec_d, frame='icrs', unit='deg')
        self.sky_main   = SkyCoord(RA_m, Dec_m, frame='icrs', unit='deg')

        # Masses and mass proxies
        self.Lx        = dbc.measurand(Lx, 'Lx', label='$L_{500,0.1-2.4}$', un='erg\\,s$^{-1}$') # X-Ray luminosity in erg/s from 0.1-2.4 keV within R500
        self.Lx_lit    = dbc.measurand(Lx_lit, 'Lx_lit', label='$L_{500,0.1-2.4}$', un='erg\\,s$^{-1}$') # X-Ray luminosity in erg/s from 0.1-2.4 keV within R500
         
        self.M200 = dbc.measurand(M200, 'M200', label='$M_{200}$', un='$M_\odot$')    # virial mass in solar masses
        self.M500 = dbc.measurand(M500, 'M500', label='$M_{500}$', un='$M_\odot$')    # M500   mass in solar masses
        self.M100 = dbc.measurand(M100, 'M100', label='$M_{100}$', un='$M_\odot$')    # M100   mass in solar masses
        self.Mvir = dbc.measurand(Mvir, 'Mvir', label='$M_{vir}$', un='$M_\odot$')    # Mvir   mass in solar masses

        """ Pre-existing electron content normalizationNorm """
        self.PreNorm   = None


        """ A general or map specific reference """
        self.reference = reference
        

        """ Maps specific properties ... in the future this could be a list of map-classes """            
        #=== Classes attached to it
        self.status   = status   # A list of possible atatuses
        self.regions  = regions  # a list of the region used to seperate relics  # One to many; regions=RelicRegion('', [], rtype=1)
        self.relics   = relics   # a list of relics One to many
        self.dinfo    = replaceNone(dinfo, DetInfo())   # detection information --> One per radio map!  ... One to one ?
        self.histo    = Histo    # Histogramm of all binned objects, anyhow a future function could provide this right out of the image(s)
        self.mockobs  = replaceNone(mockobs, MockObs(0))  # mockobs   information --> One per radio map! One to one ?
        self.compacts = []  # A list of (compact) sources that were substacted. One to many

        #== 
        self.halo = halo                                     

        #=== Image --> Library of Image arrays with Image (2D numpy array or astropy object with metainformation & labels)
        # an dictionary of 2D numpy array, attached to frequencies, mock observations and other stuff
        self.mapdic     = mapdic   # a list of paths to the [vanilla] or [vanilla, subtracted, subtracted_convolved, residuum] 
        self.spixel     = 0        # has to become filled the pixle size --> however pixelsize is also part of detinfo, ... at the end I might merge an map-class and detinfo class
#        self.center    = []       # has to become filled with two center coordiantes
        self.massproxy_mask = []   # a mask that can be used to mask the radio image, if it has the same size like the image that was used to create the mask. So ... it has many ifs
        self.ClassFlag  = int(ClassFlag)
        
        
        ## All Fluxes are in mJy
        self.flux_lit   = dbc.measurand(flux_lit, 'F_lit', un='mJy', label='$F_\mathrm{lit}$')    # Literature flux
        self.Prest100   = Prest100
        #  contaminating sources:
        self.contSour  = []         # A list of contaminating sources, The source themself could be represented by objects   
       
        
        """ filter criteria """
        self.maxcomp    = 0.17   # shape criterion
        self.updateInformation(massproxis=True)           # Information with added value:

 
    """
    functions
    """
    def set_center(self): 
        """ This function does not have any current use.
        
        input: sky coordinates
        
        At its completed state it should be used to redefine the cluster centre.
        This however is only of use if the the polary binned histogramms are adjusted to that.
        At this step is currently done on the fly, using this function is a thing of future versions """

    def maps_update(self, nparray, mapname, fitsname, dinfo=None):
        """ This functions is used to write a .fits file to disk. and to update the map dictionary
        
        Input arguments:
            The numpy array map
            The mapname
            The fitsname
            The detection information
            
        Currently, there is a detour by copying the dictionary.
        I did this to counter the effect, that the dictionaries of all galaxy clusters ended up do be the same.
        """
        if dinfo is None:
            dinfo = self.dinfo
        
        mapdic = copy.deepcopy(self.mapdic)
        folder = os.path.dirname(fitsname)
        if not os.path.exists(folder):
            os.makedirs(folder)
        fitsut.map2fits(nparray, dinfo, fitsname)   # s_pixel, center1, center2  
        mapdic[mapname] = fitsname
        self.mapdic = mapdic

    #==== importing dill and 'dilling' the objects instead of pickling can help     
    # The lambda function does make the file unpickelable
    #  self.Filter = lambda x: (x.flux > self.minflux and (x.iner_rat/(x.LAS/(x.dinfo.beam[0]/60.))) < self.maxcomp)  

    def filterRelics(self, Filter=True, maxcomp=None, regard=[1,2,3]):

        if Filter:
            minflux = self.dinfo.rms * 8 * 1000  #microjansky to millijansky
            if maxcomp is None:
                maxcomp = self.maxcomp
        else:
            maxcomp = 1e9
            minflux = -1
            
        return [relic for relic in self.relics
                if ((relic.flux() > minflux) and (relic.region.rtype.classi in regard))] #and (relic.shape_advanced().value < maxcomp)

    def add_regions(self, regions, **filterargs):
       
        if len(regions) > 0:
            self.regions = self.regions + regions
            self.updateInformation(**filterargs)

    def add_relics(self, relics, **filterargs):
        
        if len(relics) > 0:
            self.relics = self.relics + relics
            self.updateInformation(**filterargs)
#     
#    def cutout_relics(self, relics = lambda x: x.relics,  Filter=False):
#        
#        # use filter to get to now which relics you extract from the maps
#        # ...
#        Image = self.Image 
#        # apply image mask based on contours        
#        
#        return Image
       
    def updateInformation(self, massproxis=False, **filterargs):
        """ Updates the cluster information, based on the associated regions/relics
            Also includes properties, which are by definition accociated to clusters with radio relics
        """     

        # Astropy cosmology
        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.27)
        self.cosmoDL = self.cosmo.luminosity_distance(self.z.value).cgs.value # cm
        self.cosmoPS = self.cosmo.kpc_proper_per_arcmin(self.z.value).value /60  # kpc/''

        # Computes M200 from other mass proxies if needed  
        if massproxis:
            self.infer_M200_MX(self.M500, 500)
            self.infer_M200_MX(self.M100, 100)
            self.infer_M200_MX(self.Mvir, 178)
            self.infer_M200_LX()

            #starts to compute all other quantities from this value
            self.infer_LX_M200()
            self.M500 = self.M(500)
            self.R200 = self.R(200)   # a function of the mass and redshift

            """ Lambda_vir = 178 for the critical density of the universe
                The critical density is very similar to the average density
            """
            self.Mvir = self.M(178)   # virial mass in solar masses
            self.Rvir = self.R(178)   # virial radius in Mpc (comoving volume?)
          
            """Development """
            if self.M200.value == 0:
                self.M200.value = np.nan

        """ relics """
        self.sparseA = np.empty(0)
        self.sparseD = np.empty(0)
        self.sparseW = np.empty(0)
        self.moment_angle = dbc.measurand(0, 'moment_angle', un='', label='$moment_\mathrm{angle}$')
        self.flux      = dbc.measurand(0, 'F', un='mJy')
        self.flux_ps   = dbc.measurand(0, 'F_ps', un='mJy', label='$F_\mathrm{ps}$')
        self.P_rest    = dbc.measurand(0, 'P_rest', un='W/Hz', label='$P_\mathrm{rest}$')
        self.Dproj     = dbc.measurand(0, 'Dproj', un='kpc', label='$D_\mathrm{proj}$')
        self.Dproj_pix = dbc.measurand(0, 'Dproj_pix', un='kpc', label='$D_\mathrm{proj,pix}$')
        self.area      = dbc.measurand(0, 'area', un='arcmin$^2$', label='$\sum\\nolimits_{i \in \mathrm{relics}} A_{i}$')
        self.area100kpc = dbc.measurand(0, 'area100kpc', un='(100kpc)$^2$', label='$\sum\\nolimits_{i \in \mathrm{relics}} A_{i}$')
        self.largestLAS = dbc.measurand(0, 'LASmax', un='kpc', label='LAS$_\mathrm{max}$')
        if self.histo is not None:
            self.histo.hist = 0*self.histo.hist
        
        # This includes the option to filter relics due to certain criteria
        for relic in self.filterRelics(**filterargs):
            relic.asign_cluster(self)
            self.area += relic.area
            self.area100kpc += relic.area * (self.cosmoPS*60/100)**2 
            self.Dproj.value = np.divide((np.multiply(self.Dproj, self.flux) +  np.multiply(relic.Dproj,relic.flux)), (self.flux+relic.flux))
            self.Dproj_pix.value = np.divide((np.multiply(self.Dproj_pix, self.flux) + np.multiply(relic.Dproj_pix, relic.flux)), (self.flux+relic.flux)) 
            self.Dproj.set_std(max(0.1*self.Dproj, 100.)) 
            self.Dproj_pix.set_std(max(0.1*self.Dproj_pix, 100.))
            self.flux += relic.flux
            self.flux.std_add(relic.flux.std)
            self.P_rest += relic.P_rest
            self.P_rest.std_add(relic.P_rest.std)
            self.flux_ps += relic.flux_ps
            self.flux_ps.std_add(relic.flux_ps.std)
            self.largestLAS.value = max(self.largestLAS(), relic.LAS())
            #self.regions.append(relic.region) """Created BUG; update needed? 
            self.area.std += [0.2*np.sqrt(relic.area/self.dinfo.Abeam[1])]*2                 # arcmin^2
            self.area100kpc.std += [area*(self.cosmoPS*60/100)**2 for area in relic.area.std]   # in (100kpc)^2
            self.accumHisto(relic.polHist)
            self.accumPixels(relic)

        self.relics_polarDistribution()
        return self

    def accumHisto(self, Hist):
        """ Adds histograms of sparse numpy arrays """
        
        if self.histo is not None and Hist is not None:
            self.histo.hist += Hist.A

    def accumPixels(self, relic):
        self.sparseA = np.concatenate((self.sparseA, relic.sparseA), axis=0)  # self.sparseA.extend(relic.sparseA)
        self.sparseD = np.concatenate((self.sparseD, relic.sparseD), axis=0)  # self.sparseD.extend(relic.sparseD)
        self.sparseW = np.concatenate((self.sparseW, relic.sparseW), axis=0)  # self.sparseW.extend(relic.sparseW)

    def comp_LX_M200(self, M=0):
        """ Boehringer+Chon+2014 , Equ. 10 """
        alpha = 1.51
        h100  = 0.70
        h70   = 1.0
        # np.sqrt(self.cosmoPara['WM']*(1+self.z)**3 + self.cosmoPara['WV'])
        Ez = np.sqrt(self.cosmo._Om0*(1+self.z)**3 + self.cosmo._Ode0)
        return 1e44 * 0.1175*Ez**alpha/h100**(2-alpha) * (M*h70/1e14)**alpha * h70**2


    def comp_M(self, m, overdensity):
        """ works with https://github.com/joergdietrich/NFW the original mass overdensity based dependent
        & depends on  z"""
        umass = u.Quantity(m, u.solMass)
        return  massC.m200_to_mdelta(umass, massC.duffy_concentration, overdensity, args=(self.z, self.cosmo)).value 


    def M(self, overdensity=200): # the original virial mass overdensity based on z
        if self.M200() > 1e5:
            result = self.comp_M(self.M200(), overdensity)
            return dbc.measurand(result, 'M%i' % overdensity, label='$M_{%i}$' % overdensity, un='$M_\odot$')
        else:
            return dbc.measurand(0, 'M%i' % overdensity, label='$M_{%i}$' % overdensity, un='$M_\odot$')
        
    def R(self, overdensity=200):
        # inspired by http://www.star.bris.ac.uk/bjm/lectures/cluster-cosmology/3-scaling.pdf
        #M200 = 4np.pi/3 * overdensity * rho_crit *  r(200)**3
        f_a = 3.0857e21  # kpc to cm
        f_b = 1.98855e33 # solar masses g to 
        rho_crit = self.cosmo.critical_density(self.z.value) #<Quantity 9.31000324...e-30 g / cm3>
        result = np.power( self.M(200)*f_b / ( 4*np.pi/3 * overdensity * rho_crit.value) ,  1./3. ) / f_a
        return dbc.measurand(result, 'R%i' % overdensity, label='$R_{%i}$' % overdensity, un='$kpc$')
        
    def infer_M200_LX(self):

        if self.M200() == 0 and self.Lx() != 0:
            masses_checked = [10**m for m in np.arange(13.5, 16, 0.01)]
            masses_tested = []
            for m in masses_checked: 
                masses_tested.append( np.abs(np.log(self.comp_LX_M200(m)/self.Lx)) )
            m_index = np.argmin(masses_tested)
            self.M200.value = masses_checked[m_index]   # [0:2pi]
          
    def infer_M200_MX(self, Mx, x):
        """ This method is VERY slow, it takes around 2 seconds per computation """

        if self.M200() == 0 and Mx() != 0:
            masses_checked = [10**m for m in np.arange(13.0, 16, 0.03)]
            masses_tested = []
            for m in masses_checked: 
                masses_tested.append(np.abs(np.log(self.comp_M(m,x)/Mx)))
            m_index = np.argmin(masses_tested)
            self.M200.value = masses_checked[m_index]   # [0:2pi]      
                      
    # inspired by Boehringer et al 2014  (also see Nuza+2017)
    def infer_LX_M200(self):

        if self.Lx() == 0 and self.M200() != 0:
            self.Lx.val = self.comp_LX_M200(self.M200)


    def flux2power_woRedening(self):
        # To factor for converting flux [Jy?] to power (erg/s)
        return 1e-29*(4 * np.pi*np.power(self.cosmoDL*1e-2, 2))

    
    def relics_polarDistribution(self, histo=None, **kwargs):
        """ This function identifies to axis of preferred radio emission in the galaxy cluster
            In the first step the dividing line that minimizes/maximizes an regression criterium.
            In the next step the fluxes within these regions are computed.
            
            This procedure depends on gcl.histo, if flux is outside the given range (radius etc)
            it won't be put into consideration
        """

        eps = 1.e-15  # 1  Femto Jy as offset ... hopefully this doesn't change the results

        if histo is None:
            histo = self.histo
        if histo is not None:
            collabsed = np.sum(histo.hist, axis=1)   # Collapses along the distance axis
            N = collabsed.shape[0]

            # This is some simple form of regression!
            fitarray_pro = [np.sum(np.multiply(np.power(collabsed, 1.0), np.abs(np.cos(histo.ticks[0]-shift)))) for shift in histo.ticks[0]]
            fitarray_anti = fitarray_pro
            self.fitarray_pro = fitarray_pro

            # The value fitted would give the counterclockwise angle. This is why I subtract them from 2pi
            self.relic_pro_index = np.argmax(fitarray_pro)
            self.relic_pro_angle = histo.ticks[0][self.relic_pro_index]      # [0:2pi]

            def find_nearest_index(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return idx

            #self.moment_angle.value = compute_moments(self.sparseA, self.sparseD, self.sparseW)
            #self.relic_pro_index = find_nearest_index(histo.ticks[0], self.moment_angle)
            #print('____________', self.moment_angle*180/np.pi, self.relic_pro_angle*180/np.pi,
            #      min(self.moment_angle.value-self.relic_pro_angle,  self.moment_angle.value-self.relic_pro_angle)*180/np.pi)
            #print('\n\n')
            #self.relic_pro_angle = self.moment_angle

            self.relic_anti_index = np.argmin(fitarray_anti)
            self.relic_anti_angle = histo.ticks[0][self.relic_anti_index]   # [0:2pi]

            self.pro1 = np.sum(maut.polar_from_to(collabsed, [int(self.relic_pro_index-1./4.*N), int(self.relic_pro_index + 1./4.*N)])) + eps
            self.pro2 = np.sum(collabsed) - self.pro1 + 2.*eps

            self.anti1 = np.sum(maut.polar_from_to( collabsed, [int(self.relic_anti_index),int(self.relic_anti_index + 1./2.*N)])) + eps 
            self.anti2 = np.sum(collabsed) - self.anti1 + 2.*eps 
    
            # simple but inefficient way to make sure that the vectors have one and not two major directions
            if self.pro2 / self.pro1 > 1:
     
                self.relic_pro_index = int((self.relic_pro_index + len(fitarray_pro)/2) % len(fitarray_pro))
                self.relic_pro_angle = histo.ticks[0][self.relic_pro_index]      # [0:2pi]
                self.pro1, self.pro2 = self.pro2, self.pro1
    
            if self.anti2 / self.anti1 > 1:
                
                self.relic_anti_index = int((self.relic_anti_index + len(fitarray_pro)/2) % len(fitarray_pro))
                self.relic_anti_angle = histo.ticks[0][self.relic_anti_index]   # [0:2pi]
                self.anti1, self.anti2 = self.anti2, self.anti1

        else:
            self.pro1, self.pro2, self.anti1, self.anti2 = eps, eps, eps, eps
            
        # Divide flux in up- and downside emission
        self.ratio_pro = dbc.measurand(self.pro2 / self.pro1, 'ratio_pro', label='ratio$_\mathrm{pro}$',  vrange=[2e-3,1])
        self.ratio_anti = dbc.measurand(self.anti2 / self.anti1, 'ratio_ant', label='ratio$_\mathrm{anti}$', vrange=[2e-3,1])
        ratio_val = mstats.gmean([self.ratio_pro(), self.ratio_anti()])
        self.ratio_relics = dbc.measurand(ratio_val, 'ratio_relics', label='ratio$_\mathrm{dipol}$',
                                           std=[abs(ratio_val-min(self.ratio_pro(),self.ratio_anti())),
                                                abs(max(self.ratio_pro(),self.ratio_anti())-ratio_val)],
                                           vrange=[2e-3,1])
    
    
    def gettypestring(self, vec=False) :
        """ Returns a typestring for the diffusive emission components in this galaxy cluster.
        In the case of several emission types of different candidate status I have to think of an solution"""
        
        tarr = [0]*4

        for reg in self.regions:
            if reg.rtype.classi == 0 or reg.rtype.classi == -1:
                tarr[0] = 1.0-0.6*int(reg.candidate)
            if reg.rtype.classi == 1 or reg.rtype.classi == 3:
                tarr[1] = 1.0-0.6*int(reg.candidate)
            if reg.rtype.classi == 2:
                tarr[1] = 1.0-0.6*int(reg.candidate)
                tarr[2] = 1.0-0.6*int(reg.candidate)
             
        if len(self.regions) > 1 and tarr[2] == 0:
            tarr[0] = 0.5
        
        if self.halo == 'TRUE':
            tarr[3] = 1.0 
        if self.halo == 'CAND':
            tarr[3] = 0.4
                
        if vec:
            return tarr
        else:
            tstr = ''
            tstr += '\\textcolor{black!%.f}{\\textbf{(}}' % (tarr[1]*100)
            tstr += '\\textcolor{black!%.f}{$\\bullet$}'  % (tarr[3]* 50)
            tstr += '\\textcolor{black!%.f}{\\textbf{)}}' % (tarr[2]*100)
            tstr += '\\textcolor{black!%.f}{*}'           % (tarr[0]*100)
            return tstr

    
    def getstatusstring(self, zmin=0.05) :
        """ Returns if the galaxy cluster was considered and the analysisi and the statusstring  descriping either the total flux density value or
            the resson why the cluster was not included in the analysisi """
        
        if self.flux > 0 and self.z > zmin:
            return True, '$%.1f$' % self.flux.value
        else:
            """ NVSS specific """
            
            status = self.status.split(',')
            statuses = []
            
            if self.z < zmin:
                statuses.append('low z')
            
            if self.flux() == 0:
                vec = self.gettypestring(vec=True)
                if 'noMAP' not in status:
                    if vec[1] > 0 not in status:
                        statuses.append('too faint' ) 
                    else:
                        statuses.append('not gischt')
                    
            if 'CLASSI'    in status: statuses.append('not gischt')
            if 'CONF'      in status: statuses.append('contam.')
            if 'FAINT'     in status: statuses.append('too faint')
            if 'noMAP'     in status: statuses.append('low Dec')
            if 'NOCLUSTER' in status: statuses.append('no cluster')



            statuses = list(set(statuses))  # removes duplicates
            statusstring = ','.join(statuses)

            return None, statusstring
    
    def where_all(self, allwhere=None, **kwargs):
        """ Returns a numpy-where list for a two dimensional array that an be used to mask all relics in an 2D-map"""

        for relic in self.relics:  # self.filterRelics(**kwargs):
            if allwhere is not None:
                allwhere = np.concatenate((relic.pmask, allwhere), axis=1)
            else:
                allwhere = relic.pmask

        return allwhere
    
    def whereToMask(self, H, iL):
        """ Returns a numpy array for a mask and an array of a given shape """
        return [iL[0].astype(np.int64), iL[1].astype(np.int64)]
    
    def Mask_Map(self, H, normalize=1, **kwargs):
        """Input:
           H       :  A numpy array that represents the values of the evenly sampled map
           Normalize: Array or Float for Normalization 
           
           Out: A masked and normalized numpy array weighted by the 'nornalization' factor
           """
        iL = self.where_all(**kwargs)
        mask = self.whereToMask(H, iL)
 
        Hnew = np.asarray(H)*np.nan
        Hnew[mask] = H[mask]/normalize

        return Hnew


    def stoch_drop(self, seed_dropout):
        """ stochasticaly drops cluster objects based on  discovery_prop_cluster()"""
        if seed_dropout is None:
            return (len(self.relics)>0) and (len(self.filterRelics()) > 0)
        return self.discovery_prop_cluster() > seed_dropout.uniform(0, 1)
    
    def discovery_prop_cluster(self):
        """ Assigns a probalistic value (0-1) that the given relic hosting cluster would be detected in any relic survey 
            In the current implementation, this is a super simple approach that assumes that
            the most prominent object determines the discovery probability.
            
            In the truth characteristic structures like double relics might increase the discovery probability.
            This is why if there is more than 2 objects above a certain value,
            we double the discovery probability ... in the future
            """

        return max(surut.discovery_prop(self.relics))
        
    def __str__(self):
        """ I don't see any reason to keep this return value. Is this for any specific filtercluster statement?"""
        return True 
    
        
# Define a new object: galaxycluster_added for the CW, in the future might as well marge bith classes
class Galaxycluster_simulation(Galaxycluster):
    """
    Simple class for representing a galaxy in a simulation. In a future implementation one could incorporate the used cosmology in this class,
    alternatively one could consider merging this with the  'mockobs' class.
    """

    def __init__(self, name, ids, z=0.05, mass_gas=0, Lx=0, R100=0, Prest_Vol=0, PrestVazza=0, ratio=0, **kwargs):
        """
        decription
        """
        super(Galaxycluster_simulation, self).__init__(name, 0., 0., z, Lx=Lx, **kwargs)
        #CBclass.Galaxycluster.__init__(self, name, 0, 0, z, Lx=Lx, dinfo=dinfo, mockobs=mockObs, Histo=Histo)

        self.ids = ids

        self.Prest_vol = dbc.measurand(Prest_Vol, 'Prest_vol', label='$P_\\mathrm{rest,vol}$', un="W/Hz")   # total radio emission within virial radius (of projected quantities), could also be a future functionality for an efficiency of 1
        self.PrestVazza = PrestVazza

        """ these measures are particular used for the CW simulations """
        self.R100_vaz = dbc.measurand(R100, 'R100_vaz', label='$R_\mathrm{100,vaz}$', un='Mpc')    # R200 in [Mpc]
        self.interratio = ratio   # log10 (Interpolation Ratio)

    def __str__(self):
        return 'Galaxycluster_simulation: Nothing specified in def__str__(self)'


class DetInfo:
    """ This class is used to define the observational parameters of an created map, albeit one survey/ each galaxy cluster could have regions with differing detaction informations
    """
    def __init__(self, name='', beam=[1,1,0], spixel=1, rms=0, limit=0, nucen=1.4, center=None, pcenter=[0,0], survey = 'NVSS', telescope='VLA-D'):

        self.FWHM2sigma = 1/2.355
        self.name    = name     # name/id of the survey , e.g. NVSS, MSSS, TGSS
    #    self.survey  = survey  # survey  as a string, e.g. NVSS, MSSS, TGSS --> outdated and now in the survey class
        self.beam    = beam     # beam size FWHM [major [arcsec], minor[arcsec], theta]
    
        self.spixel  = spixel   # pixel size in arcsec
        self.rms     = rms      # rms noise in Jy/beam
        self.limit   = limit    # detlimit  in Jy/beam 
        self.nucen   = nucen    # observed frequency in GHz
     
        self.center    = center
        self.pcenter   = pcenter
        self.telescope = telescope
        
        self.update_Abeam()

    def update_Abeam(self):
        self.beam_pix= [self.beam[0]/self.spixel, self.beam[1]/self.spixel, self.beam[2]]                       # beam size [major [pix], minor[pix], theta]
        self.Abeam   = [1.133*self.beam_pix[0]*self.beam_pix[1], 1.133*self.beam[0]*self.beam[1]/3600] 

    def convolve_map(self, npmap):
        """Convolves a map and converts the  to numpy quantities.
        The sum of the array is not preserved, because of the way the map represent flux per beam"""                      
        return self.Abeam[0]*ndi.gaussian_filter(npmap, (self.FWHM2sigma*self.beam_pix[0], self.FWHM2sigma*self.beam_pix[1]))  ## gaussian convolution
                
               
    def __str__(self):
        # return("The cluster %12s at dec=%6.2f and rec = %6.2f has an z of %5.3f" % (self.name, self.dec, self.rec, self.z))
        return 'Nothing specified in def__str__(self) for <detinfo>'


class BFieldModel:
    """
    Simple class for representing the parametrisation of the magnetic field
  
    """             

    def __init__(self, B0=0, kappa=0, compress=0.85, B_para='dens'):
    
        # Bott variables used in the formula
        # B = (B0/muG) * (rho/1e-4) ** kappa
        self.B0       = B0        # B0
        self.kappa    = kappa     # kappa
        self.compress = compress  # factor of magnetic field enhancement in case of compression. 0..1
        self.B_para   = B_para       # 'dens' or 'press'. Governes to what the B-field parameter is referring to

    def __str__(self):
        return 'B-field model: B0=%6.3e muG, kappa=%6.3f' % (self.B0, self.kappa)


class RModel(BFieldModel):
    """
    Class that inherits from 'BFieldModel'  and adds parameters to descripe the model parameters for DSA
    Also has the parameters for a population of preexisting electrons and further, misc parameters.
    """

    def __init__(self, id, effList = [1], simu=True, pre=False, **kwargs):

        BFieldModel.__init__(self, **kwargs)

        """ Metaparameter """
        self.pre     = pre
        self.effList = effList
        self.simu    = simu
        self.id      = id

class PreModel_Hoeft(RModel):
    """
    Class that inherits from 'RModel'  and adds parameters to descripe a population of preexisting electrons and further, misc parameters
    """
    def __init__(self, id, t0=0.5, t1=7, n0=1e-6, n1=1e-1, ratio=1e-5, **kwargs):

        """ Model of preexisting electrons, see internal .pdf description"""
        RModel.__init__(self, id, pre=True, **kwargs)
        self.t0 = t0   # Minimal time since reacceleration in Gyr
        self.t1 = t1   # Maximal time since reacceleration in Gyr
        self.n0 = n0   # Electron number density of accretion shocks
        self.n1 = n1   # Electron number density of 'core'
        self.ratio = ratio  # Initial normalisation PRE and thermal at shockfront

class PreModel_Gelszinnis(RModel):
    """
    Class that inherits from 'RModel'  and adds parameters to descripe a population of preexisting electrons and further, misc parameters
    """

    def __init__(self, id, p0=1e-4, p_sigma=1, sigmoid_0=1e-4, sigmoid_width=1, **kwargs ):

        """ Model of preexisting electrons, see internal .pdf description"""
        RModel.__init__(self, id, pre=True, **kwargs)
        self.p0            = p0             # p_rat
        self.p_sigma       = p_sigma        # scatter sigmoid in lognormal coordinates
        self.sigmoid_0     = sigmoid_0      # density  (particles/cm-3)value for which the sigmoid value becomes 0
        self.sigmoid_width = sigmoid_width  # order of magnitude over which this plasma content changes
        self.PREexpo       = 0.09           # The exponent for the eff(Mach) modification, latter this should be replaced by the average gamma of the PREs
        
class MockObs:
    """
    Simple class for representing some parameters of the mock observation
    Development: could also be an descendant if DetInfo
    """

    def __init__(self, ii, theta=0, phi=0, psi=0, proj='', snap=-1, z_snap=0, clid=-1, hsize=6000, binmax=1000, headerc='Mpc', snapfolder='', xrayfolder=''):

        self.id     = int(ii)
        self.theta  = theta
        self.psi    = psi
        self.phi    = phi

        self.proj   = proj

        self.hsize  = hsize   # in kpc, gives the Diameter of the Cutout
        self.binmax = binmax  # in pixel

        self.z_snap = z_snap
        self.snap   = int(snap)
        self.clid   = int(clid)

        self.snapfolder = snapfolder  # Folder for the original snapshots
        self.xrayfolder = xrayfolder  # Folder for the original snapshots
        self.headerc    = headerc

    def __str__(self):
        return 'MockObs class object with id: %i has  Snap %i & ClID %i, angles=(%6.3f,%6.3f,%6.3f) in radiants' % (self.id, self.snap, self.clid, self.theta, self.psi, self.phi)


   
class DetRegion(object):

  """
  Simple class for representing a relic region in a FITS image.
  --> In future it just just be a DetectionRegion with an assigned object class
  """             

  def __init__(self, name, cnt, cnt_WCS=[]):
    
     
     """ DEVELOPMENT BEGIN """
     self.name       = name
     self.detections = []
     self.cnt        = cnt    # In list of an list ... (I know!] cnt coordinates
     self.cnt_WCS    = cnt_WCS
    
#     def accumulate(classi=['relics']):   
#         for detections in self.detections('class'=relics):
#             """accumulate properties, like for galaxy clusters"""
#             LAS  =1
#             flux =1
     
      
     """ DEVELOPMENT END """


class RelicRegion(DetRegion):  # Define a new object: RelicRegion

  """
  Simple class for representing a relic region in a FITS image.
  --> In future it just just be a DetectionRegion with an assigned object class
  """             

  def __init__(self, name, cnt, dinfo=DetInfo(), rtype=-1, alpha=None, alpha_err=False, alphaFLAG=False, candidate=False, **kwargs):

    super(RelicRegion, self).__init__(name, cnt,  **kwargs)   

    self.rtype     = RelicTypes(rtype)
    self.candidate = candidate
    self.alphaFLAG = alphaFLAG
    
    if alphaFLAG:
       self.alpha = -np.abs(float(alpha))  #We assume the input alpha to be negative
    else:
       self.alpha = float("NaN")
                             
                             
    if alpha_err:
      try:
       self.alpha_err = float(alpha_err)
      except: 
       self.alpha_err = 0.3
    else: 
       self.alpha_err = 0.3
       
       
    """DEVELOPMENT
    #'' Maps pecific properties ... in the future this could be a list of map-classes ''#       
    #=== Classes attached to it
    self.status     = status   # A list of possible atatuses
    self.relics     = relics  # a list of relics One to many  
    self.dinfo      = dinfo   # detection information --> One per radio map!  ... One to one ? 
    self.histo      = Histo    # Histogramm of all binned objects, anyhow a future function could provide this right out of the image(s)
    self.mockobs    = mockobs # mockobs   information --> One per radio map! One to one ?
    self.compacts   = compacts # A list of (compact) sources that should be / were substacted            
                                           
    #=== Image --> Library of Image arrays with Image (2D numpy array or astropy object with metainformation & labels)
    self.Image      = Image    # an library of 2D numpy array, attached to frequencies and mock observations  --> outdated? as it should be saved in fitsfiles?
    self.fitsname   = fitsname # a list of paths to the [vanilla] or [vanilla, subtracted, subtracted_convolved, residuum]
    self.xname      = xname    # a list of paths to x-ray images
    self.spixel     = 0        # has to become filled the pixle size --> however pixelsize is also part of detinfo, ... at the end I might merge an map-class and detinfo class
#        self.center    = []       # has to become filled with two center coordiantes
    self.massproxy_mask = []   # a mask that can be used to mask the radio image, if it has the same size like the image that was used to create the mask. So ... it has many ifs
    self.ClassFlag  = int(ClassFlag)  
    
    ## All Fluxes are in mJy
    self.flux_lit   =  dbc.measurand(flux_lit, 'F_lit' , un = 'mJy' , label='$F_\mathrm{lit}$')    # Literature flux
    self.Prest100   = Prest100
    self.PrestVazza = PrestVazza
 
    #  contaminating sources:
    self.contSour  = []         # A list of contaminating sources, The source themself could be represented by objects   
   
    
    #'' filter criteria ''#
    self.maxcomp    = 0.17   # shape criterion
    self.minflux    = minflux    #.6    # mJy
    self.updateInformation(massproxis=True)           # Information with added value:
            

            
    DEVELOPMENT END"""  
            
  def __str__(self):
    return("The RelicRegion %s of relic type=%3i with an alpha of %5.3f" % (self.name, self.rtype, self.alpha))
    
class RelicTypes:  

    def __init__(self, classi = -1):

        # Add the plotting linestyles to them
        # Make each type an own class that inherits from diffuse emission
        self.classi         = classi
        self.classes        = ['unspecified', 'Phoenix', 'gischt', 'double gischt', 'gischtlet'][self.classi+1]
        self.classessimple  = ['-', 'PhX', 'RL', 'dRL', 'SL'][self.classi+1]
        self.alpha          = [-1.25, -1.8, -1.25, -1.25, -1.25][self.classi+1]
    
    def __str__(self):
        # return("The cluster %12s at dec=%6.2f and rec = %6.2f has an z of %5.3f" % (self.name, self.dec, self.rec, self.z))
        return ( 'rtype(%i -> %s)' % (self.classi, self.classessimple) )
        #return("The cluster %12s at dec=%11s and ra = %9s has an z of %5.3f" % (self.name, self.dec, self.ra, self.z))
      
class ApperanceAtFrequ(patches.Ellipse):
    """ REDUNTANT in current implementation"""
    def __init__(self, flux=1, frequ=1.4, major=0, minor=0, posangle=0):
          self.flux   = flux  # [mJy at 1.4 GHz]
          self.frequ  = frequ
          self._width = major
          self._height = minor
          self._angle = posangle
    
    
class RadioSource(ephem.FixedBody):
    """ REDUNTANT in current implementation"""
    def __init__(self, ra=0, dec=0, flux=1, frequ=1.4, ha= 0, major=0, minor=0, posangle=0):  #ra='00:00:00', dec='00:00:00'
        ephem.FixedBody.__init__(self)
        self._ra  = ra
        self._dec = dec
        self.ha = ha
        self.apper = [ApperanceAtFrequ(flux, frequ, major, minor, posangle=0)]
    
    def AddObs(self, apperance):
        self.apper.append(apperance)
      

class Relic:
    """
    A class for representing an observed relic (radio relic, phoenix, etc.) in our visible universe.
    To initialize all values it needs a galaxy cluster to be assigned
    ... region and dinfo might be outdated
    ... ! the angle computation might be needed to  updated !
    ... make sure that all possible output quantities that are not arrays, lists ore class objects are measurands
    ... inconsistencies exist for the alpha value and its usage
    ... might add astropy functionalities --> including using astropy quantities (quantity + unt)
    ... In future this might become a child class of the RadioSource class
    """

    def __init__(self, region, dinfo, RA, Dec, LAS, area, GCl=False, F=0, F_ps=0, F_lit=0, LLS_lit=0, alpha=None,
                 alpha_err=0, alphaFLAG=False, theta_elong=False, cov=False, cnt=[], Mach=None, Dens=0, Dproj_pix=0,
                 polHist=None, sparseD=None, sparseW = None, sparseA=None, pmask=None):
        """
        parameters
        """
        self.GCl    = None      # galaxycluster
        self.name   = ''
        self.region = region    # RelicRegion that defines an subarray of the detection region, as well as labels and  (as a proxy) spectral information
        self.dinfo  = dinfo     # Detection Info Class
        self.cnt    = cnt       # Contours in the detection image given as pixel coordinates (?)
        
        #===  Hist
        self.polHist = polHist  # a sparse numpy array

        ##===  Averaged: Position and coordiantes
        self.RA  = dbc.measurand(RA, 'RA', un='deg')
        self.Dec = dbc.measurand(Dec, 'Dec', un='deg')
          
        self.Dproj_pix = dbc.measurand(Dproj_pix, 'Dproj_pix', label='$D_\\mathrm{proj,pix}$', un="kpc", std=max(0.1*Dproj_pix, 100.))
                                     
        ##=== Geometry and Morphology
        self.LAS      = dbc.measurand(LAS, 'LAS', un="'", std=0.4)   # '
        self.LAS_beam = dbc.measurand(LAS/(dinfo.beam[0]/60), 'LAS_beam', label='LAS$_\\mathrm{beam}$', un="$beam$", std=0.4/(dinfo.beam[0]/60))   # '
        self.area     = dbc.measurand(area, 'area', label='$A_\\mathrm{relic}$', un="arcmin^2", std=0.2*np.sqrt(area/dinfo.Abeam[1]))
        
        #=== All Fluxes are in mJy  
        self.flux_ps = dbc.measurand(F_ps, 'F_ps', label='$F_\\mathrm{ps}$', un='mJy', std=0.10*F_ps)  #shoud be a function!
        flux_err     = np.sqrt((dinfo.rms*1e3*self.area()/dinfo.Abeam[1])**2 + (0.05*self.flux_ps())**2 + (0.10*F)**2)
        self.flux    = dbc.measurand(F, 'F', un='mJy', std=flux_err)
                    
                    
        #=== Literature values
        self.flux_lit = dbc.measurand(F_lit, 'F_lit', label='$F_\\mathrm{lit}$', un='mJy', std=0.10*F_ps)  #F_lit  #Flux at 1.4 GHz
        self.LLS_lit  = dbc.measurand(LLS_lit, 'LLS_lit', label='LLS', un='kpc', std=max(0.1*LLS_lit, 100.))   #sLLS_lit*1000
          
        #=== spectral index and k correction
        self.alpha     = dbc.measurand(alpha, 'alpha', '$\\alpha_\mathrm{int}$', un=None, std=alpha_err)
        self.alphaFLAG = alphaFLAG
        self.spec_cor = None
        
        # Flux error needs area
        
        self.filling     = self.area()/(np.pi/4*self.LAS()**2)
        self.filling_err = self.filling*np.sqrt((self.area.std[0]/self.area())**2 + (2*self.LAS.std[0]/self.LAS())**2)
             
        self.cov         = cov        #2x2 covariance matrix
        self.theta_elong = theta_elong
        
        """ Shape information """
        eigvals, eigvecs = np.linalg.eigh(cov)
        self.eigvecs    = np.sort(eigvecs, axis=1)[::-1]   #not yet working perfectly! I think there is a bug inside, you could test by multiplying eigvecs with egvalues to see if the covariance matrix is reproduced!
        self.eigvals    = np.sort(eigvals)[::-1]
        self.iner_rat   = dbc.measurand(self.eigvals[1]/self.eigvals[0], 'iner_rat', label='$v_\\mathrm{PC2}/v_\\mathrm{PC1}$', un=None, vmax=1)
        self.iner_rat_n = dbc.measurand(self.eigvals[1]/self.eigvals[0], 'iner_rat_n', label='$v_\\mathrm{PC1}/v_\\mathrm{PC2}$', un=None, vmin=1)
        self.ecc        = np.sqrt(1-self.iner_rat())  # BUG shoud work with out ()

        """ Arrays: Position and coordiantes """
        self.pmask   = pmask     # The ordering is due to the pixels in the fits array! I.e. pixel positions can be reconstructed
        self.sparseD = sparseD   # distance
        self.sparseW = sparseW   # flux weight
        self.sparseA = sparseA   # angle
         
        """ pmask is used to fasten the binning of following quantities"""                
        self.wT = None  # a sparse numpy array; flux weighted Temperature
        self.wRho = None  # a sparse numpy array; flux weighted electron upstream density
        self.wMach = None  # a sparse numpy array; flux weighted Mach number
        self.wArea = None  # a sparse numpy array; flux weighted Mach number
        self.wAlpha = None  # a sparse numpy array; flux weighted Mach number
        self.wB = None  # a sparse numpy array; flux weighted Mach number
        self.wPre = None  # a sparse numpy array; flux weighted Mach number
        self.wDoG_rel = None
#        self.wRho_up     = None
#        self.wRho_down   = None
#        self.wT_up       = None
#        self.wT_down     = None

        self.averages_quantities()
        self.infer_Mach(Mach)  
        
        if isinstance(GCl, Galaxycluster): 
            self.asign_cluster(GCl)
        else:
            'This relic has no galaxy cluster assigned!'

    """
    FUNCTIONS
    """

    # Average mach number, from https://arxiv.org/pdf/1611.01273.pdf (2016) ; page 13, equ. 5; though they have an inconsistentdefinitionof alpha(assume alpha>0 in equation)
    # Note that: This equation is tailored for alphas that are defined to be negative  for relics.
    #             2alpha+4/ 2alpha
    #              -2alpha-inj+3/(-2alpha_inj+1)
    #               -2alpha+4/(-2alpha+2)
    #               -alpha+2/-alpha+1

    """ See Colafrancesco+2017 equ. 5 & 6 """
    def infer_Mach(self, Mach=None):
        """ Derive the Mach number proxy if the spectral index is known """

        #          print (self.alpha())
        if Mach is not None:
            self.Mach = dbc.measurand(Mach, '$\overline{M}$', un=None)
        else:
            if self.alpha.value is None:
                Mach = np.nan  #np.sqrt( (-alpha+1) / (-alpha-1 ) )  #=
            else:
                alpha = min(self.alpha.value, -1.03)
                Mach = np.sqrt((-alpha+1) / (-alpha-1))
            self.Mach = dbc.measurand(Mach, '$\overline{M}$', un=None)

      
    def create_Histo(self, GCl, normtype='R200'):
        """ Creates a PolarHistogram if the Histogram() type is set """

        if GCl.histo is not None:

            Nexp = GCl.histo.norm.Nexp
            normtype = GCl.histo.norm.ntype

            if normtype == 'R200':
                norm = GCl.R200() * (1500/GCl.R200())**(Nexp-1)
            if normtype == 'Rvir':
                norm = GCl.Rvir() * (1500/GCl.Rvir())**(Nexp-1)
            if normtype == 'Dproj':
                norm = 1e3

            self.polHist = sparse.csr_matrix(np.histogram2d(self.sparseA, self.sparseD / norm, bins=GCl.histo.bins,
                                                            normed=False, weights=self.sparseW)[0])

    def shape_advanced(self):
        return dbc.measurand(self.iner_rat/ (self.LAS / (self.dinfo.beam[0]/60.)), 'shape_advanced',
                             label='shape$_\mathrm{PC_1}$', un=None)
                  

    def averages_quantities(self):
        """ flux weighted properties """
        if self.wMach is not None:
            self.T       = dbc.measurand(np.sum(self.wT)/np.sum(self.sparseW), 'T_av', label='$\overline{T}$', un='keV')
            self.Mach    = dbc.measurand(np.sum(self.wMach)/np.sum(self.sparseW), 'M_av', label='$\overline{M}$', un=None)
#            self.Rho     = dbc.measurand(np.sum(self.wRho)/np.sum(self.sparseW), 'rho_av', label='$\overline{\\rho_\mathrm{averaged}$', un='cm$^{-3}$') #label='$\overline{\rho}_\mathrm{down}$',
            self.Area_av = dbc.measurand(np.sum(self.wArea)/np.sum(self.sparseW)*1e6, 'area_av', label='$\overline{\mathrm{area}}$', un='kpc$^2$')
            self.Area    = dbc.measurand(np.sum(self.wArea)/np.sum(self.sparseW)*len(self.sparseW), 'area', un='Mpc$^2$')
            self.B       = dbc.measurand(np.sum(self.wB)/np.sum(self.sparseW), 'B', un='$\mu G$')
            self.fpre    = dbc.measurand(np.sum(self.wPre)/np.sum(self.sparseW), '$f_\mathrm{Pre}$', un=None)
            self.DoG_rel = dbc.measurand(np.sum(self.wDoG_rel), 'DoG_\mathrm{rel}')
            self.Rho_up  = dbc.measurand( np.sum(self.wRho_up)/np.sum(self.sparseW), 'rho_up', label='$\overline{\\rho}_\mathrm{up}$', un='cm$^{-3}$')
#                self.Rho_down = dbc.measurand( np.sum(self.wRho_down)/np.sum(self.sparseW)          , 'rho_dow', label='$\overline{\\rho}_\mathrm{down}$', un = 'cm$^{-3}$' )
#                self.T_up     = dbc.measurand( np.sum(self.wT_up)/np.sum(self.sparseW)              , 'T_up'   , label='$\overline{T_\\mathrm{up}}$', un = 'keV' ) 
#                self.T_down   = dbc.measurand( np.sum(self.wT_down)/np.sum(self.sparseW)            , 'T_dow'  , label='$\overline{T_\\mathrm{down}}}$', un = 'keV' ) 

        if self.wAlpha is not None:
            self.alpha = dbc.measurand(-np.sum(self.wAlpha)/np.sum(self.sparseW), 'alpha', label='$\\alpha$', un=None)
#                self.infer_Mach()

    def asign_cluster(self, GCl):
        self.GCl = GCl
        self.create_Histo(GCl)
        if self.region.name == '':
            self.name = GCl.name
        else:
            self.name = GCl.name + '_' + self.region.name

        self.LLS        = dbc.measurand(self.LAS * 60 * GCl.cosmoPS, 'LLS', un="kpc", std=[std*60*GCl.cosmoPS for std in self.LAS.std])
        self.area100kpc = dbc.measurand(self.area*(GCl.cosmoPS*60/100)**2, 'area100kpc', un="(100kpc)$^2$", std=[std*(GCl.cosmoPS*60/100)**2 for std in self.area.std])

        # All radio luminosities are in 10^24 W/Hz, factor 1e-29 -->  1e-3 Jy/mJy  *  (1e-26 W/(Hz m^2)) / Jy
        # All radio luminosities are in W/Hz
        # See Nuza et all. 2012
        self.bw_cor = 1./(1.+GCl.z.value)     # Bandwidth quenching not covered by luminousity distance
        P           = self.bw_cor * self.flux.value*GCl.flux2power_woRedening()
        self.P      = dbc.measurand(P, 'P', label='$P_\\mathrm{obs}$', un="W/Hz",
                                      std=P*np.sqrt( (self.flux.std[0]/self.flux.value)**2 + ( self.region.alpha_err*np.log(1+GCl.z.value)*(1+GCl.z.value)**self.alpha_z() )**2) )  # include z error
        self.P_rest = dbc.measurand( self.P.value*self.speccor(GCl), 'P_rest', label='$P_\\mathrm{rest}$', un = "W/Hz",
                                      std=[std*self.speccor(GCl) for std in self.P.std])   #!!! generall alpha;
        self.P_lit  = dbc.measurand( self.flux_lit()*1e-29*(4* np.pi*np.power(GCl.cosmoDL * 1e-2, 2)), 'Plit', label='$P_\\mathrm{lit}$', un = "W/Hz")

        self.vec_GCl   = [(GCl.RA.value-self.RA())*np.cos(self.Dec()*np.pi/180), self.Dec.value-GCl.Dec.value]       # corrected for spherical coordiante system
        self.theta_GCl = (np.angle(self.vec_GCl[0]+self.vec_GCl[1]*1j, deg=True)+360) % 360  # in deg;  The complex number form  allows for np.angle. Theta increases counterclockwise

        self.Dproj     = dbc.measurand( np.linalg.norm(self.vec_GCl)*3600*GCl.cosmoPS, 'Dproj', label='$D_\\mathrm{proj}$', un='kpc')
        self.Dproj_rel = dbc.measurand( np.linalg.norm(self.vec_GCl)*3600*GCl.cosmoPS, 'Dproj', label='$D_\\mathrm{proj}$', un='kpc')

        if self.Dproj > 1e4:
            self.Dproj.value = 0

        if self.theta_elong:  # [ Ellipsecontours, Moment of Inertia, Ellipse moments ]
            self.theta_elong = self.theta_elong                                                  # in deg; vector measured     np.angle
            self.theta_rel   = dbc.measurand(maut.MinAngle_quarter(self.theta_elong, self.theta_GCl), 'theta_rel', label='$\\theta_\mathrm{relic-GCl}$', un='$^\circ$', vrange=[0,90])   # in deg; vector measured     np.angle
        else:
            self.theta_elong = False
            self.theta_rel   = False

    def alpha_z(self, default=True):
        if default or self.alpha is None:
            alpha = self.region.rtype.alpha
        else:
            alpha = self.alpha()
        return alpha

    def speccor(self, GCl, default=True):
        speccor = 1./np.power(1.+GCl.z(), self.alpha_z(default) )
        self.spec_cor = speccor
        return speccor
        
    def __str__(self):
        return 'relic.name: %s' % (self.name)