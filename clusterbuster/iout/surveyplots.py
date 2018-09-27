#!/usr/bin/env python

''' 
Created on 2017
@author: jakobg

This .py provides the functionalities of higher-level data products from the pickled relic catalogues, including:
- Creating pandas data tables (implement saving them as .csv or o.ods files)
- Creating .pandas scatter matrixes
- Creating .fits or .png images of the simulated objects and their galaxy clusters
- Creating ...
'''

from __future__ import division,print_function




import copy
import os
import warnings
import aplpy # check https://aplpy.readthedocs.io/en/v0.9.9/_generated/aplpy.aplpy.FITSFigure.html  for comands

import clusterbuster.surveyclasses as cbclass
import clusterbuster.dbclasses     as dbc 
import clusterbuster.iout.misc     as iom
import clusterbuster.maput         as maput
import matplotlib.pyplot           as plt
import matplotlib.colors           as colors
import math                        as math
import numpy                       as np

from   matplotlib      import cm
from matplotlib.ticker import NullFormatter
from scipy import stats
import matplotlib.patches as patches



import surveysim.music2.mockobsxray as Xray
from PyPDF2 import PdfFileMerger

def plot_RelicEmission_polar(surveys, compsurvey=None, single=True, modeltext=True, additive=False, 
                             aligned=False, cbar=True, addinfo=False, mirrored=False,minrel=1,eff=None, 
                             zborder  = 0.05, plottype = 'flux', title="Polar binned radio relic flux", 
                             dpi=None, Histo=dbc.Histogram2D()):
    ''' posible inprovements: http://stackoverflow.com/questions/22562364/circular-histogram-for-python 
    minrel : minimal number of relics to be consided for the histogramm

    if surveys is a list of surveys this function will plot averaged quantities

    '''
    development = True # Some features in development

    
    import collections  # to assert if we deal with an list of surveys or a single object
    if not isinstance(surveys, collections.Iterable):   surveys = [surveys]
    

    if surveys[0].mainhist is not None:
        Histo = surveys[0].mainhist

    
    
    ''' Plotting and normalization of combined histogramm'''
    expPlot  = 0.45   # 0.25
    cmap     = cm.viridis # ; cm.afmhot none
    ztype    = '>'
    
    yticks   = [0.4,0.8,1.2]
    add      = 1./2.  # This angle is added anticlockwise! We want to turn by 90 degree anticlockwise
    addangle = int(aligned)*np.pi*add  # Additional Rotation               

    halfHists = []
    radials   = []
    stats     = []
      
    if compsurvey  is not None: 
         _ , comprad, comppol, _ , _= compsurvey.polar(zborder = zborder,ztype = ztype)
         comprad = compsurvey.polar(zborder = zborder,ztype = ztype)[1][0]    #comprad[1]
         deviations= [] 
          
    for survey in surveys:

        survey.set_binning(Histo)
        
        nowfolder = '%s/Relics_polar/' % (survey.outfolder) 
        iom.check_mkdir(nowfolder)

        buckedfolder = os.path.abspath(os.path.join(survey.outfolder, '..', 'bucket'))
    #    buckedfolder = '../%s' % (survey.outfolder) 
        iom.check_mkdir(buckedfolder)
        effs = survey.Rmodel.effList
        
        for eff in effs:
#            print (eff) #DEBUGGING
#            for ii,GCl in enumerate(survey.FilterCluster(minrel=minrel,eff=eff,zborder=zborder,ztype='>')):
#                
#
#                ''' Deriving the Histogramm should be a functionality of the survey or the relic cluster, so this should become outdated '''
#                GCl.updateInformation(eff=eff, Filter=True)  
#                if GCl.histo is not None and np.sum(GCl.histo.hist) != 0:
#                    inner  = Histo.bins[1][0:-1]
#                    outer  = Histo.bins[1][1::]
#                    angle  = Histo.ticks[0]
#                    angles, innerZ = np.meshgrid(angle, inner, sparse=True)
#                    angles, outerZ = np.meshgrid(angle, outer, sparse=True)
#                    AreaHist  = 2*np.pi*(2*np.pi)/len(angle)*(outerZ**2-innerZ**2)   
#                    shiftHist  = np.roll(Histo.hist.T, -int(aligned*(GCl.relic_pro_index)), axis=1) / AreaHist**(survey.expA)  #+1e-10        
#                                        
#                    # Plots the single clusters
#                    if single:
#                        fig, ax    = plt.subplots(figsize=(14,14),subplot_kw=dict(projection='polar'), dpi=dpi)  #,subplot_kw=dict(projection='polar')  
#                        ax.pcolormesh(Histo.bins[0], Histo.bins[1],  shiftHist, cmap=cmap)
#                        ax.set_theta_offset(addangle)
#                        ax.annotate("", xy=(int(not aligned)*(GCl.relic_pro_angle), 1.5), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))
#                        from matplotlib import transforms as mtransforms
#                        ax.arrow(0, 0, 0, 0.5, linewidth=3, width=0.005, transform=mtransforms.Affine2D().translate( int(not aligned)*(GCl.relic_pro_angle), 0) + ax.transData)
#            
#                            
#                        ax.text(0.01, 1.05, '%s' %(GCl.name.replace('_',' ')), fontsize=20, transform=ax.transAxes)
#                        if addinfo:
#                            ax.text(0.3, 0.9, 'Summed relic flux: %.2e Jy' %(np.sum(Histo.hist)), fontsize=20, transform=ax.transAxes, color='w')
#                            ax.text(0.3, 0.87, 'Ratio pro: %.2e  anti: %.2e' %(GCl.ratio_pro(), GCl.ratio_anti()), fontsize=20, transform=ax.transAxes, color='w') 
#                        ax.set_yticks(yticks) 
#                        ax.tick_params(axis='x'                , labelsize=25)
#                        ax.tick_params(axis='y', colors='white', labelsize=25)
#                        ax.grid(True)
#                        if title is not None: ax.set_title(title, va='bottom')
#                        for ftype in ['pdf','png']:
#                            plt.savefig('%s%s-polar.%s' % (nowfolder,GCl.name, ftype))
#
#                        fig.clf()     
                        
    
            
            if additive:   
                _, (radial,radial_tickz), halfHist, stat, mesh = survey.polar(eff=eff, aligned=True, minrel=minrel,  mirrored=mirrored, mode=plottype, zborder = zborder,ztype = ztype)
              
                if halfHist is not None:
                    fig, ax    = plt.subplots(figsize=(14,14),subplot_kw=dict(projection='polar'), dpi=dpi)  #,subplot_kw=dict(projection='polar')
                    if not mirrored:
                        meshed = ax.pcolormesh(Histo.bins[0][0:int(len(Histo.bins[0])/2)+1], Histo.bins[1],  halfHist, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot)  ) #, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-1.0, vmax=1.0),cmap='RdBu_r'                         
                        if compsurvey: ax.pcolormesh(Histo.bins[0][int(len(Histo.bins[0])/2):], Histo.bins[1],  comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot)  )                 
                    else:   
                        meshed = ax.pcolormesh(Histo.bins[0][int(len(Histo.bins[0])/2):], Histo.bins[1],  halfHist, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot)  )
                        if compsurvey: ax.pcolormesh(Histo.bins[0][0:int(len(Histo.bins[0])/2)+1], Histo.bins[1],  comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot)  )                  
                    ax.set_theta_offset(addangle)
                    if addinfo: ax.text(0.3, 0.9,  'Summed relic flux: %.2e Jy (all cluster)' %(stat), fontsize=20, transform=ax.transAxes, color='w')
                    ax.set_yticks(yticks)
                    ax.tick_params(axis='x',                 labelsize=25)
                    ax.tick_params(axis='y', colors='white', labelsize=25)
                    ax.grid(True)
                    if cbar:  fig.colorbar(meshed)
                    if title is not None: ax.set_title(title, va='bottom')
                    for ftype in ['pdf','png']: #,'jpg'
                        plt.savefig('%s/%s%.2f-polar.%s' % (nowfolder,survey.name, np.log10(eff), ftype))
                #        HistoDD('Output/%s/PostProcessing/AllTogether' % (survey.name), Histo, survey.name, logging=False)
                    fig.clf()   
    

                halfHists.append(halfHist)
                radials.append(radial)
                stats.append(stat)
                if compsurvey is not None: deviations.append(np.sum(np.abs(radial-comprad)))
#
#            if development:
#                fig, ax    = plt.subplots(figsize=(14,14), dpi=dpi)  #,subplot_kw=dict(projection='polar') 
#                ax.pcolormesh(np.asarray(range(35))*3.5/35, np.asarray(range(46))*1.5/46, mesh, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot)  )  # Histo.ticks[0]/2
#                plt.xlabel('$R_{200}/\mathrm{Mpc}$')
#                plt.ylabel('$D_\mathrm{proj}/R_{200}$')
#                plt.savefig('%s/butterfly-%s%.2f.%s' % (   nowfolder,survey.name, np.log10(eff), ftype))
#                
                
                # plot ratio of relics flux, 
                # plot average/median pro relic distance
                # plot sigma pro rleic distnace
                # plot skew
                
                
    
    ''' Colorsheme from seaborn '''
#    cmap = ListedColormap(sns.color_palette('deep'))
    if len(radials) > 0:

        
        ''' Radial plot '''
        fig, (ax1)    = plt.subplots(1, 1, figsize=(7,4), dpi=dpi)  
        ax1.set_xlabel('Distance [$R_{200}$] along pro-relic axis')
        ax1.set_ylabel('Weighted signal W')
    
        if compsurvey is not None: 
            ax1.plot(radial_tickz, comprad, alpha=0.6, color='blue')  
            
        for radial in radials:
            ax1.plot(radial_tickz, radial, alpha=np.sqrt(1/len(radials)),c='grey')   # color=cmap(0.0)
        
        
        # Made up confidence intervals (I'm too lazy to do the math...); math comes later
#        variance = radials
        ''' Development '''
#        var        = np.var(np.asarray(radials), axis=1)
#        mid_bound  = np.mean(np.asarray(radials))
#        high_bound = mid_bound+var
#        low_bound  = mid_bound-var
#                    
#        print(high_bound.shape, var),            
#        ax1.fill_between(radial_tickz, high_bound, low_bound, color='gray', alpha=0.4)
#        ax1.plot(radial_tickz, radial, color='gray', alpha=0.9)   # color=cmap(0.0)


        ax1.set_xticks(yticks + [0] + [-y for y in yticks]) 
        ax1.set_xlim(-Histo.bins[1][-1],Histo.bins[1][-1])
        dist = 0.08
        
        
        ''' Textlabels '''
        if modeltext and survey.Rmodel.simu:
            mod = survey.Rmodel 
            kwargs = {'verticalalignment':'bottom', 'horizontalalignment':'right', 'transform':ax1.transAxes, 'color':'black', 'fontsize':15, 'alpha':0.8}
            ax1.text(0.40, 0.90       , '$\log_{10}(eff) =%+.2f$,' % ( np.log10(eff)),   **kwargs)
            ax1.text(0.40, 0.90-1*dist, '$\log_{10}(B_0) =%+.2f$,' % ( np.log10(mod.B0)),**kwargs)
            ax1.text(0.40, 0.90-2*dist, '$\kappa       = %+0.2f$ ' % ( mod.kappa),       **kwargs)
            if compsurvey  is not None: 
                ax1.text(0.40, 0.90-3*dist, '$\Delta\\,\\mathrm{signal}  = %0.2f$ '% ( np.average(deviations)),       **kwargs)
            if isinstance(mod, cbclass.PreModel_Hoeft):
                ax1.text(0.40, 0.90-4*dist, '$t_{1;2}  = %0.3f\,;\,%0.3f$ '% ( mod.t0, mod.t1) ,       **kwargs)
                ax1.text(0.40, 0.90-5*dist, 'ratio$\\mathrm{_{pre}}  = %0.2f$ '% ( mod.ratio ),       **kwargs)                
    
            
            if survey.Rmodel.pre:
                ''' NOT implemented yet '''
#                print( p0, pre, p_sigma, sigmoid_0, sigmoid_width )        
    #            ax2.set_yscale('log')
        for ftype in ['pdf','png']: #,'jpg'
            plt.savefig('%s/%s%.2f-sumprofile.%s' % (   nowfolder,survey.name, np.log10(eff), ftype))
            plt.savefig('%s/%s%.2f-sumprofile.%s' % (buckedfolder,survey.name, np.log10(eff), ftype))

        
        # Statistics, contribution
        fig, (ax1)    = plt.subplots(1, 1, figsize=(7,4), dpi=dpi)  
        
        for stat in stats:
            statistic = np.divide(stat,np.sum(stat))
            ax1.hist(statistic, color='b', alpha=1/len(stats), bins='auto')  # arguments are passed to np.histogram
        for ftype in ['pdf','png']: #,'jpg'
            plt.savefig('%s/%s%.2f-sumsstats.%s' % (   nowfolder,survey.name, np.log10(eff), ftype))
            plt.savefig('%s/%s%.2f-sumsstats.%s' % (buckedfolder,survey.name, np.log10(eff), ftype))
            
        
        ''' Polar plot'''
        halfHist = np.sum(halfHists,axis=0)/len(halfHists)
        fig, ax    = plt.subplots(figsize=(14,14),subplot_kw=dict(projection='polar'), dpi=dpi)  #,subplot_kw=dict(projection='polar')
        nearmax    = np.partition(halfHist.flatten(), -2)[-2]
        if not mirrored:
            meshed = ax.pcolormesh(Histo.bins[0][0:int(len(Histo.bins[0])/2)+1], Histo.bins[1],  halfHist, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=nearmax  ) #, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-1.0, vmax=1.0),cmap='RdBu_r'                         
            if compsurvey is not None: 
                     ax.pcolormesh(Histo.bins[0][int(len(Histo.bins[0])/2):], Histo.bins[1],   comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=nearmax  )                 
        else:   
            meshed = ax.pcolormesh(Histo.bins[0][int(len(Histo.bins[0])/2):], Histo.bins[1],  halfHist, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=nearmax  )
            if compsurvey is not None: 
                     ax.pcolormesh(Histo.bins[0][0:int(len(Histo.bins[0])/2)+1], Histo.bins[1], comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=nearmax  )
        ax.set_theta_offset(addangle)
        if addinfo: ax.text(0.3, 0.9,  'Summed relic flux: %.2e Jy (all cluster)' %(stat), fontsize=20, transform=ax.transAxes, color='w')
        ax.set_yticks(yticks)
        ax.tick_params(axis='x',                 labelsize=25)
        ax.tick_params(axis='y', colors='white', labelsize=25)
        ax.grid(True)
        if cbar:  fig.colorbar(meshed)
        if title is not None: ax.set_title(title, va='bottom')
        for ftype in ['pdf','png']: #,'jpg'
            plt.savefig('%s/%s%.2f-polar-sums.%s' % (   nowfolder,survey.name, np.log10(eff), ftype))
            plt.savefig('%s/%s%.2f-polar-sums.%s' % (buckedfolder,survey.name, np.log10(eff), ftype))
    #        HistoDD('Output/%s/PostProcessing/AllTogether' % (survey.name), Histo, survey.name, logging=False)
    
        '''DEBUGGING'''
        fig.clf()   

            
    return 0
        
def plot_Clusters(survey, dynamicscale=False, subtracted=True, relicregions=False, DS9regions=False, diamF=2.6, 
                  colorbar=False, beam=True, shapes=False, recenter=True, infolabel = False, sectors=False,
                  xray=False, highres=False, show_rot=False, vectors = False, label_sheme='balanced', 
                  maxdist=1700, filterargs = {'zborder':0, 'ztype':'>', 'minimumLAS':4, 'GClflux':20, 'index':None}):
     
    pdfs  = []
    eff   = survey.Rmodel.effList[0]
#    cmap  = 'afmhot'   #'BuPu'   #
#    norm  = MPLcolors.PowerNorm(gamma=2) #CW uses a gamma of 2
    #cmap =plt.cm.get_cmap('RdYlBu')
    laargs = {'color':'#BBBBBB'}      # line arguments
    ciargs = {'color':'#BBBBBB'}      # arguments for the circle/centre area
    baargs = {'color':'#BBBBBB'}      # argument for the scale bar
    if label_sheme =='dark':
        laargs.update({'color':'black'})
        ciargs.update({'color':'#111111'})
        baargs.update({'color':'#111111'})
    elif label_sheme =='bright':    
        laargs.update({'color':'w'})
        ciargs.update({'color':'w'})
        baargs.update({'color':'w'})


    for GCl in survey.FilterCluster(**filterargs):
    
        GCl = copy.deepcopy(GCl) # Because else changes will influence all galaxy clusters, you knwo class referencing in python, i.e.
#        print(GCl.filterRelics(eff))
#        print(filterargs)
#        self.Rmodel.effList[0]
#        print('!!!', [relic.flux() for relic in GCl.filterRelics(eff)])
#        print(GCl.largestLAS(),GCl.flux(), GCl.filterRelics(eff), GCl.filterRelics(), eff)
#        print('!!!', GCl.relics)  
#        return 0

 
        if xray:  
            '''X-Ray''' 
            if 'Brems' not in GCl.mapdic:     
                savefolder = survey.outfolder
                Xray.Run_MockObs_XRay(GCl, savefolder, verbose=False)  
                ''' Is the mapdic updated ? '''
        else:
            GCl.mapdic['Brems'] = GCl.mapdic['Diffuse']
        
        #=== Sets some plotting parameters
        kpc = GCl.cosmoPS*3600  # kiloparsec per degree
         
        R200 = GCl.R200()  
        if R200 <= 0:  
            print( 'For GCl %s no mass proxy and hence no virial radius is known. For plottting issues we set R200=1600 kpc.' % (GCl.name) ) 
            R200 = 1600
        radius  =       R200/kpc
        diam    = diamF*R200/kpc 
 
        import seaborn as sns; sns.set(style="white", color_codes=True) #ticks
        f = aplpy.FITSFigure(GCl.mapdic['Brems']) #dimensions=[0, 1],, slices=[10, 10], , slices=[10, 10]
#        f.add_grid()
#        f.grid.hide()
#        f.grid.set_xspacing(0.2)
#        f.grid.set_color('blue') 
#        f.grid.show()

#        f.remove_grid()

        if recenter:
            print('Recentering', GCl.name, GCl.RA(),GCl.Dec(),diam)
            f.recenter(GCl.RA(),GCl.Dec(),width=diam,height=diam) # radius is also possible!

        if survey.Rmodel.simu:
#            f.tick_labels.set_xformat("dd:mm:ss")
#            f.tick_labels.set_yformat("dd:mm:ss")  
            f.axis_labels.set_xtext('Coordinate 1')
            f.axis_labels.set_ytext('Coordinate 2')
            f.tick_labels.hide()
        else:
            f.tick_labels.set_xformat("dd:mm:ss")
            f.tick_labels.set_yformat("dd:mm:ss") 
            f.axis_labels.hide()



        # The basic image
        if dynamicscale:
            vmax    = np.max(f._data)
            levels  = [GCl.dinfo.limit*l for l in [ survey.m_cnt**n for n in np.arange(0,16)  ]]  #0,16
        else:
#            vmid    = 3e-7
#            exponent= 1.7 
            print('plot_Clusters fixed scale might be brocken!')
            vmax    = survey.emi_max
            levels  = survey.cnt_levels
           

        vmin    = 0.6 * GCl.dinfo.rms #0.25
        vmid    = -2  #0.00006  * GCl.dinfo.rms
        exponent= np.log(vmax/vmin)
         
#        print(levels, survey.cnt_levels)
#        levels   = [l*1e-3 for l in levels]


        if not xray:
            
            for relic in GCl.filterRelics(eff=eff):
                pixelcnt = np.transpose(np.squeeze(relic.cnt))
#                print( pixelcnt.shape )
                wcscnts  =  f.pixel2world(pixelcnt[0,:],pixelcnt[1,:])
#                print( type(wcscnts), len(wcscnts) )
                wcscnts =  np.asarray([ (x,y)  for x,y in zip(wcscnts[0],wcscnts[0]) ]).T
#                print( wcscnts.shape )
                f.show_polygons([wcscnts], lw=2, color = 'white') # , lw=1, color = 'white' , alpha=1.0
#                print( len(relic.cnt) )
            addargs = {'vmid':vmid,'vmin':vmin,'vmax':vmax,'stretch':'log','exponent':exponent}
            
            ''' It seems like you can only have one interactive contours '''

#            f.show_colorscale(vmin=1e9, vmax=1e11,  stretch='linear', cmap='afmhot') #gist_heat 
            f.show_colorscale(vmid=vmid, vmin=vmin, vmax=vmax,  stretch='log', exponent=exponent, cmap='afmhot') #gist_heat      
#            f.show_contour(GCl.mapdic['Diffuse'], linewidth=0.15, overlap = True, levels=[l for l in levels if l<vmax*survey.m_cnt], cmap='afmhot',filled=True, alpha=0.49, extend='max', **addargs) #
            print(levels, survey.cnt_levels)
            f.show_contour(GCl.mapdic['Diffuse'], linewidth=0.15, overlap = True, levels=levels, colors='green',filled=False)
        else:
             if 'MUSIC' in survey.name or 'Threehundret' or 'ShockTest' in survey.name: #MUSIC-2
                vmin_xr    = 2.5
                vmax_xr    = 9.7 #6.2
                vmid_xr    = -1.5
             else:
                vmin_xr    = -2
                vmax_xr    =  5.
                vmid_xr    = -4.5
            
             exponent= np.log(max(vmax/vmin,1.0001))
            
             f.show_colorscale(vmid = vmid_xr, vmin=vmin_xr, vmax=vmax_xr,  stretch='log', exponent=exponent, cmap='afmhot') #gist_heat      
#             f.show_colorscale(cmap='afmhot') #gist_heat  
#             return 0
#             f.show_contour(GCl.xname[-1], colors='grey', linewidth=0.5,  levels=levels, overlap = True)
             #development
#             ''' X-ray '''
#        levels = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7]
#        levels = [l-3.5 for l in levels]   
                         
             if highres:
                 key      = "Raw"
                 key_comp = "CompModell" 
                 levels = [levels[0]/8, levels[0]/4,levels[0]/2] + levels
             else:
                 key      = "Diffuse"
                 key_comp = "Subtrated"
                 

#             addargs = {'vmid':vmid,'vmin':vmin,'vmax':vmax,'stretch':'log','exponent':exponent}
             if key_comp in GCl.mapdic:
                 f.show_contour(GCl.mapdic[key], linewidth=0.15, overlap = True, levels=levels, colors='green',filled=False)
             if key_comp in GCl.mapdic and subtracted: 
                 f.show_contour(GCl.mapdic[key_comp], linewidth=0.15, overlap = True, levels=levels, colors='red', filled=False )  

             
        
        if shapes:
            for relic in GCl.filterRelics(maxcomp=100) : 
                vlen    =  (np.sqrt(relic.iner_rat())*relic.LLS+0.05*R200)/kpc
                f.show_arrows(relic.RA(), relic.Dec(), -relic.eigvecs[1][0]*vlen, relic.eigvecs[1][1]*vlen) #-np.cos(relic.Dec*np.pi/180)*
                f.add_label(relic.RA-relic.eigvecs[1][0]*vlen, relic.Dec+relic.eigvecs[1][1]*vlen, '$s = %.2f$' % (relic.iner_rat()), size='x-large', **laargs)

        # The Jakobs circle OR (virial) radius
        f.show_circles(GCl.RA(), GCl.Dec(), radius, linestyle='--', **ciargs)
        if sectors:
            GCl.relics_polarDistribution(histo=survey.mainhist)
            P1   = [GCl.RA() - np.cos(GCl.relic_pro_angle)*radius*1.05/np.cos(np.radians(GCl.Dec())),GCl.Dec() + np.sin(GCl.relic_pro_angle)*radius*1.0]
            P2   = [GCl.RA() + np.cos(GCl.relic_pro_angle)*radius*1.05/np.cos(np.radians(GCl.Dec())),GCl.Dec() - np.sin(GCl.relic_pro_angle)*radius*1.0]
            
            P1b   = [GCl.RA() - np.cos(GCl.relic_anti_angle)*radius*1.0/np.cos(np.radians(GCl.Dec())),GCl.Dec() + np.sin(GCl.relic_anti_angle)*radius*1.0]
            P2b   = [GCl.RA() + np.cos(GCl.relic_anti_angle)*radius*1.0/np.cos(np.radians(GCl.Dec())),GCl.Dec() - np.sin(GCl.relic_anti_angle)*radius*1.0]
            f.show_lines( [ np.array(zip(P1,P2)) ]  , color='w', lw=2., linestyle=':')
            f.show_lines( [ np.array(zip(P1b,P2b)) ], color='r', lw=2., linestyle=':') 
            
            if GCl.ratio_relics() > GCl.ratio_relics.vrange[0]:  # Plot if multiple relic
                f.add_label(P1[0],P1[1],  'ratio= %.1e' % (GCl.ratio_relics()), size='x-large', **ciargs)
        
        # Workaround ... in future it would be better to take the image information from the image and read the contours directly
        
        try:
            _, center, spixel = maput.FITS2numpy( GCl.mapdic['Raw'] )
        except:
            _, center, spixel = 0,(0,0),7.5

        if relicregions:
            #contours, contourinfos = iom.readDS9regions('Regions/RR_%s.reg'% (GCl.name), spixel, center[0], center[1], pixelcoords=False)  
            styles  = ['--',':','-','-','-']
            f.show_polygons( [np.transpose(np.squeeze(region.cnt_WCS)) for region in GCl.regions], lw=2, linestyle=styles[GCl.relics[0].region.rtype.classi+1], **laargs) # , alpha=1.0, facecolor='orange'
        
        if DS9regions: # Load a regions file into APLpy plot
            f.show_regions('Regions/RR_%s.reg'  % (GCl.name))
           
        f.add_scalebar(1)
        f.scalebar.show(1000./kpc, linestyle='solid', linewidth=3., alpha=0.7, **baargs)  
        f.scalebar.set_corner('bottom left')
        f.scalebar.set_label('1 Mpc')
        f.scalebar.set_font(size='large', weight='medium', stretch='normal', family='sans-serif',  style='normal', variant='normal')

        if beam:
            f.add_beam()
            f.beam.show()
#            f.beam.set_major(GCl.dinfo.beam[0] * u.arcsecond)
#            f.beam.set_minor(GCl.dinfo.beam[1] * u.arcsecond)
#            f.beam.set_angle(GCl.dinfo.beam[2])  # degrees
            f.beam.set(facecolor='#BBBBBB',alpha=0.8,edgecolor='black')

        if show_rot:     
            f.add_label(0.97, 0.12, '$\\theta=%.3f$' % (GCl.mockobs.theta),  relative=True, style='oblique', size='large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
            f.add_label(0.97, 0.08, '$\\phi  =%.3f$' % (GCl.mockobs.phi)  ,  relative=True, style='oblique', size='large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
            f.add_label(0.97, 0.04, '$\\psi  =%.3f$' % (GCl.mockobs.psi)  ,  relative=True, style='oblique', size='large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
                  
        if infolabel:
            f.add_label(0.97, 0.95, '  %s'     % (GCl.name.replace('_', ' ')), relative=True, style='oblique', size='xx-large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
            f.add_label(0.97, 0.90, '$z=%.2f$' % (GCl.z()),  relative=True, style='oblique', size='xx-large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
        else:
            f.add_label(0.97, 0.95, '$z=%.2f$' % (GCl.z()),  relative=True, style='oblique', size='xx-large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
            f.add_label(0.97, 0.90, '$z_\mathrm{snap}=%.2f$' % (GCl.mockobs.z_snap),  relative=True, style='oblique', size='xx-large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
            
            
        if colorbar:
            f.add_colorbar()
            f.colorbar.show()
#            f.colorbar.set_axis_label_text('%s flux density [Jy/beam]' % (survey.name_short))
            f.colorbar.set_axis_label_text('$\log_{10}(P_\\mathrm{Brems,bol}$ in restframe) [arbitrary unit]')
            

        '''DEVELOPMENT'''
        if vectors:

           pdata = GCl.mapdic['MassSpeed'] #fits with magnitude of signal (use where not enough signal)
           adata = GCl.mapdic['MassAngle'] #fits with angle of signal     (use nan, where no vector should be)
           f.show_vectors(pdata, adata, step=15,scale=1e-2,alpha=0.2, color='blue',lw=2) # , mutation_scale=4 ,ls='-.-', 0.3, head_width=5

#           x  = GCl.mapdic['x'] + 3.0 dsdsd+ GCl.RA  #fits with magnitude of signal (use where not enough signal)
#           y  = GCl.mapdic['y'] + GCl.Dec #fits with angle of signal     (use nan, where no vector should be)
#           dx = GCl.mapdic['dx']  #fits with magnitude of signal (use where not enough signal)
#           dy = GCl.mapdic['dy']  #fits with angle of signal     (use nan, where no vector should be)
#           f.show_arrows(x, y, dx, dy, step=15,scale=1e-2,alpha=0.2, color='blue',lw=2) # , mutation_scale=4 ,ls='-.-', 0.3, head_width=5
        '''DEVELOPMENT END'''
#


        nowfolder = '%s/Images/' % (survey.outfolder)
        iom.check_mkdir(nowfolder)
        savefile = '%s/%s-%s%s' % (nowfolder,survey.name, GCl.name,'HR'*highres)
        f.save('%s.png' % savefile, dpi=400)    
        f.save('%s.pdf' % savefile)
#        circular_cutout(f, savefile)
        pdfs.append('%s.pdf' % savefile)
        f.close()
        plt.close("all")
        print(len(pdfs))
        
    merger = PdfFileMerger()
    
    for pdf in pdfs:
        merger.append(pdf)
    
    merger.write('%s/%s-All%s.pdf' % (survey.outfolder,survey.name,'HR'*highres))
        


def circular_cutout(f, savefile):
    
    """
    Demo of image that's been clipped by a circular patch.
    """
    
    
#    image = plt.imread(image_file)
#    
    fig, ax = plt.subplots()
#    im = ax.imshow(image)
    patch = patches.Circle((260, 200), radius=200, transform=ax.transData)
    f._figure.set_clip_path(patch)
    
    ax.axis('off')
    f.save(savefile+'_circular.png')

#============== Reads relic information out of an ds9 region file
def readDS9relics(regfile, spixel, center, pixref, Test=False): 

   contours   , contourinfos = iom.readDS9regions(regfile, spixel, center, pixref)     
   contoursWCS, _            = iom.readDS9regions(regfile, spixel, center, pixref, pixelcoords=False)    

   rinfo =   []
   for ii,info in enumerate(contourinfos):
     #info       =  info.split(', ')
    
     try:
       info       =  info.split(', ')             # get infolist
       info[2]    =  info[2].replace('alpha ','') #remove this alpha string
       split      =  info[2].split(' ')           #split alpha into substrings
       if len(split) > 1:
          alpha      = split[0]
          alpha_err  = split[1]
       else:
          alpha      = split[0]
          alpha_err  = 0

       reg        =  cbclass.RelicRegion(name = info[0], cnt=[contours[ii]], cnt_WCS=[contoursWCS[ii]], rtype=int(info[1]), alpha=alpha, alphaFLAG=('false'  not in alpha.lower()), alpha_err=alpha_err)
     except:
       reg        =  cbclass.RelicRegion(name =  '',  cnt=[contours[ii]], cnt_WCS=[contoursWCS[ii]], rtype=-1, alphaFLAG=False )  
       
     if ('test'  not in info[0].lower() ) or not Test:
        rinfo.append( reg )
     
   return rinfo

   

def PlotDistribution_FluxRatioLAS(location, ClList, RList):
#=== Test the deviation of the  literatur evalue and measured fluxes of the brightest objects in a cluster
   
    # use latex for font rendering
    import matplotlib as mpl
    mpl.rcParams['text.usetex'] = True

    plt.clf()
    fig = plt.figure(figsize=(8, 4.7), dpi=200)
    ax  = plt.subplot(111) 
  
    # brightrelics = ....
    # filter brightest object ... get only its flux ...
    # compare it with flux minus other fluxes'
  
  
    extrapolL = [] #n fact its: ['A1443','A1682','MCS_J1149.5+2223','MCS_J2243.3-0935','ZwCl_2341+0000']
    for o in ClList:
       if  o.name in extrapolL:
           print( o.name  )
  
    x_list   = [o.largestLAS                 for o in ClList if o.flux_lit >0 and o.name not in extrapolL]
    y_list   = [o.flux/o.flux_lit            for o in ClList if o.flux_lit >0 and o.name not in extrapolL] #[np.log10(o.flux/o.flux_lit)  for o in ClList if o.flux_lit >0]
    y_err    = [o.flux_err/o.flux_lit        for o in ClList if o.flux_lit >0 and o.name not in extrapolL]
    area     = [np.power(o.flux,0.35)*4e0    for o in ClList if o.flux_lit >0 and o.name not in extrapolL]
  
    x_list2  = [o.largestLAS                 for o in ClList if o.flux_lit >0 and o.name     in extrapolL]
    y_list2  = [o.flux/o.flux_lit            for o in ClList if o.flux_lit >0 and o.name     in extrapolL] #[np.log10(o.flux/o.flux_lit)  for o in ClList if o.flux_lit >0]
    y_err2   = [o.flux_err/o.flux_lit        for o in ClList if o.flux_lit >0 and o.name     in extrapolL]
    area2    = [np.power(o.flux,0.35)*4e0     for o in ClList if o.flux_lit >0 and o.name    in extrapolL]
  
  
    farnsw_x         = [0,     0.5,   1.0,    1.5,    2.0,    2.5,    3.0,    3.5,    4.0,    4.5,    5.0,    5.5,    6.0,    6.5,    7.0,    7.5,    8.0,    8.5,    9.0,    9.5,   10.0,   10.5,   11.0,   11.5,   12.0,   12.5,   13.0,   13.5,   14.0,  14.5,   15.0,   15.5,   16.0,   16.5,   17.0,   17.5,   18.0,   18.5,   19.0,    19.5,   20.0]
    farnsw_dec74_pix = [32,     32,    32,     24,     22,     24,     27,     31,     32,     38,     42,     49,     58,     69,     89,    116,    154,    206,    273,    341,    418,    494,    570,    641,    706,    766,    820,    868,    912,    952,   986,   1016,   1042,   1066,   1085,  1101, 1114, 1127, 1136, 1143, 1148]
    farnsw_dec18_pix = [32,     32,    32,     31,     30,     30,     30,     30,     35,     39,     45,     50,     70,     96,    132,    178,    232,    293,    363,    435,    508,    581,    652,    723,    776,    832,    880,    922,    960,    993,  1022,   1047,   1069,   1088,   1106,  1118, 1130, 1140, 1147, 1153, 1157]

    farnsw_dec18 = [(1171.-y)/1140. for y in farnsw_dec18_pix]
    farnsw_dec74 = [(1171.-y)/1140. for y in farnsw_dec74_pix]
  
    colors   = ['b','g']  
    alphas   = [1-(err/val) for (val,err) in zip(y_list, y_err)]
    
    print( alphas )
    rgba_colors = np.zeros((len(alphas),4))
    # for red the first column needs to be one
    rgba_colors[:,2]  = 1.0
    # the fourth column needs to be your alphas
    rgba_colors[:, 3] = alphas

    print( alphas )
    [ax.scatter( x_list , y_list,    s=area ,  alpha=    0.60,  color=  colors[0], zorder=2)        ]#,     
    #l_i_err  = [ax.errorbar(x_list , y_list, yerr=y_err,  alpha=    0.60,  ecolor=colors[0], zorder=2, marker='+') ]#, 
    
                #ax.scatter(x_list2, y_list2, s=area2, alpha=0.35,color=colors[1])]    
    [ax.plot(farnsw_x, [y for y in farnsw_dec18], alpha=0.7, c='grey', zorder=1), 
                ax.plot(farnsw_x, [y for y in farnsw_dec74], alpha=0.7, c='grey', zorder=1)] #[ax.plot(farnsw_x, [np.log10(y) for y in farnsw_y], alpha=0.7)]   
    #leg1     = plt.legend( l_i, ['$S_\\mathrm{1.4,\\,lit}$','$S_\\mathrm{1.4,\\,lit}$ (extrapolated)'], ncol=2, frameon=True, fontsize=11,handlelength=1, loc = 3, borderpad = 0.4, handletextpad=0.2, framealpha=0.60, scatterpoints = 1) # title='log$_{10}(P_{1.4})\\,\mathrm{[W/Hz]}$'
    ax.fill_between(farnsw_x, [y for y in farnsw_dec18], [y for y in farnsw_dec74], color='grey', alpha='0.3',zorder=1)

    spixel   = 7.5
    size     = np.linspace((60+350)/spixel, (1400++350)/spixel, 30, endpoint=True)   # size array in arcsec FWHM
    H        = np.loadtxt('../clusterbuster/AdvancedImaging/ImageInnerPos-exp.out') 
    
  
    powers   = [3,30,300]  #[10**a for a in [0.5,1.0,1.5,2.0] ]
    legl     = [np.power(po,0.38)*4.5e-0                                               for po  in powers]
    l_iii    = [ax.scatter([],[], s=leg, edgecolors='none', alpha=0.6,color=colors[0]) for leg in legl  ]
    labels   = ['%i' %(powe) for powe in powers]    
    leg2     = plt.legend(l_iii, labels, ncol=4, frameon=False, fontsize=9, handlelength=1, loc = 1, borderpad = 0.4, handletextpad=0.2, framealpha=0.70, title='$S_\\mathrm{1.4,\\,NVSS}\\,\mathrm{[mJy]}$', scatterpoints = 1)
    #plt.gca().add_artist(leg1)
  
  


    ax.set_xlim(0,20.0)
    ax.set_ylim(0,1.4) #ax.set_ylim(-1,0.5)
    ax.set_xticks(np.arange(min(ax.get_xlim()), max(ax.get_xlim())+0.5, 3.0))
    ax.set_xlabel('$\\mathrm{LAS}\,[\\mathrm{arcmin}]$') #('$\\mathrm{max(LAS)}\,[\\mathrm{arcmin}]$')
    ax.set_ylabel('$S_\\mathrm{1.4,\\,NVSS} / S_\\mathrm{1.4,\\,lit}$') #ax.set_ylabel('$\\mathrm{log_{10}}(F_\\nu / F_{\\nu, lit})\\,\\mathrm{[mJy]}$')
    
    ax.set_aspect(1.0/ax.get_data_ratio())
     
     
    # 4 linien
    # beam         0.75
    # 1Mpc z=0.10  8.98
    # 1Mpc z=0.06 14.28 
    # 1Mpc z=0.05 16.16
    # nominal Largest imagable angular scale by VLA configuration 16.94
    scales    = [0.75,                                         8.98,                           16.16,                  16.94]
    textl     = ['$\\theta_\\mathrm{FWHM}$','$\\theta_\\mathrm{z=0.10}$', '$\\theta_\\mathrm{z=0.05}$', '$\\mathrm{\\theta_{VLA,D}}$']  #['$\\theta_\\mathrm{FWHM}$','$\\theta_\\mathrm{1\\,Mpc,z=0.10}$', '$\\theta_\\mathrm{1\\,Mpcz=0.06}$', '$\\mathrm{max(\\theta_{VLA,D})}$'] 
    color     = ['black','b','b', 'black']
    height    = [ 0.14, 0.2, 0.2, 0.14]
    mod       = ( (0,0), (0,0), (0,0), (0,0) )             #( (0.6,-0.06), (0,0), (0, -0), (0.6,-0.06) )  #( (0.8,-0.0), (-0.8,0), (-0.8, -0), (0.8,-0.0) )
    
    #c     = [plt.cm.rainbow( (np.log10(np.log10(m))-ax1.get_ylim()[0])/abs(ax1.get_ylim()[1]-ax1.get_ylim()[0]) ) for m in mach]
    for ii,m in enumerate(scales):
        ax.plot(  [m]*2 , [ax.get_ylim()[0], height[ii] ], '-', c=color[ii], lw=1.8, linestyle=':', alpha=0.7 ) 
        ax.text(  m-0.4+mod[ii][0], height[ii]+0.01+mod[ii][1], textl[ii], fontsize=10, color='black', alpha=0.7)
        
        
     
    weirdcases = [o.name     for o in ClList if (o.flux_lit >0 and np.log10(o.flux/o.flux_lit) > np.log10(1.3))]
    print( 'weirdcases:', weirdcases )
 
    #ax.text(0.33, 43.7, 'Clusters not included/sampled\n by few particles', fontsize=14)
    #ax.text(0.05, 45.2, 'Cluster res. $\\gg$\n NVSS resolution', fontsize=14)
    #ax.text(0.56, 45.15, '$V_\\mathrm{sim} \\gg$ $V_\\mathrm{obs}$', fontsize=14)


  
  
    #===
    plt.savefig(location+'.pdf') 
    plt.close(fig)
  

def stats_lineregress(name, data_x, data_y, verbose = False):
    
    if len([np.log10(x) for x in data_x]) > 0:
        slope, intercept, r_value, p_value, std_err = stats.linregress( [np.log10(x) for x in data_x], [np.log10(y) for y in data_y])
        if verbose:
            print('Survey %s regression --> slope  %.3e, intercept %.3e, r_value %.3e, p_value %.3e, std_err %.3e'  % (name,  slope, intercept, r_value, p_value, std_err) )
            print("r-squared:", r_value**2)
        #--> get the full sample, for this puprose create an array in z' that you can latter flatten!
        
        return slope, intercept
    
    return None, None


def create_Mass_redshift2( SurveySamples, zrange,colors, markers = np.asarray(['.','s']), effi=[], log=[False,False], logplot=[True,True], lockedaxis=False):

    """
    Takes as inputs;
    Surveysample - a list of galaxyclusters
    parA         - a lambda function to a measurand attribute of a class
    parB         - a lambda function to a measurand attribute of a class
    **args       - optional parameter
    Returns a scatter plot (a fig object) for further modification: linear or logarithmic, in the variables (you have to specify) and saves it under the variable specific name
    
    inspired by : http://www.astrobetter.com/blog/2014/02/10/visualization-fun-with-python-2d-histogram-with-1d-histograms-on-axes/  
             and: http://matplotlib.org/examples/pylab_examples/scatter_hist.html
    """
    
    # (r.iner_rat/(r.LAS/(r.dinfo.beam[0]/60.)))
    if len(effi)==0: effi = SurveySamples[0].Rmodel.effList
                 
    parA = lambda x: x.z
    parB = lambda x: x.M200
    
    cm = plt.cm.get_cmap('RdYlBu')
    z_symbols = ['o','>']
    z_snaps   = [0.0,0.5]  
            
    for eff in effi: 
        
        #======= Old style (begin)
#        fig = plt.figure(facecolor="white", figsize=(5.6, 4.2)) # aspect='equal'
#        axScatter = fig.add_subplot(1,1,1)
#        axScatter.tick_params(which='both', direction='in')


        '''======= New style (begin)'''
        
        nullfmt = NullFormatter()         # no labels
        
        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left + width + 0.02
        
        rect_scatter = [left, bottom, width, height]
        rect_histx   = [left, bottom_h, width, 0.2]
        rect_histy   = [left_h, bottom, 0.2, height]
        
        # start with a rectangular Figure
        fig = plt.figure(1, figsize=(6, 6))
        
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        # ticks inside
        axScatter.tick_params(which='both', direction='in')
        axHistx.tick_params(which='both', direction='in')
        axHisty.tick_params(which='both', direction='in')
        
        
        #show grid
        axScatter.grid()
        axHistx.grid()
        axHisty.grid()
        
        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        '''======= New style (end)'''
#        axScatter.set_xlim( (0,0.3) )
#        ax.set_ylim( (2e-2,1) )
#        ax.set_autoscale_on(False)
    
        if logplot[0]: axScatter.set_xscale('log')
        if logplot[1]: axScatter.set_yscale('log')
    

        plotl  = []
        zlist  = []
        Data_x = []
        Data_y = []
        SurveyHistkwargs = []
    #   sc = plt.scatter(xy, xy, c=z, vmin=0, vmax=20, s=35, cmap=cm)
    #   cbar = plt.colorbar(sc)
    #   cbar.solids.set(alpha=1)
    
        for Survey in SurveySamples:
            
            print('create_LAS_shape::Survey.name:', Survey.name)

            
            if Survey.Rmodel.simu: 
                feff = eff
            else:
                feff = 1
            
            if Survey.name == 'NVSS':
                Survey.scatterkwargs.update( color='black', ecolor='black', alpha=0.4, fmt=None)
                Survey.histkwargs.update( color='black', alpha=0.4)
                del Survey.scatterkwargs['color']
                del Survey.histkwargs['color']
            else:
                Survey.scatterkwargs.update( color='r', ecolor='r', fmt=None)
                Survey.histkwargs.update( color='r')
                del Survey.scatterkwargs['color']
                del Survey.histkwargs['color']
            del Survey.scatterkwargs['fmt']
                # make the right colors

            
            for ii,z_snap in enumerate([z_snaps[0]]):     
                    data_x = []; data_y=[]
#                    for List_z in [GCl for GCl in Survey.FilterCluster(minrel=0, eff=eff) if GCl.mockobs.z_snap==z_snap] :
                    List_z = [GCl for GCl in Survey.FilterCluster(minrel=0, eff=feff)] 
                    if len(List_z) > 0:  
                            data_x = [parA(meas).value for  meas in List_z]
                            data_y = [parB(meas).value for  meas in List_z]
                            proxyA,proxyB = parA(List_z[0]), parB(List_z[0]) 
            
                    Data_x.append(data_x)   
                    Data_y.append(data_y)   
                    SurveyHistkwargs.append(Survey.histkwargs)
              
        
            if Survey.Rmodel.simu: 
                for ii,z_snap in enumerate([z_snaps[0]]):  
                    data_x = []; data_y=[]
#                    for List_z in [GCl for GCl in Survey.FilterCluster(minrel=0, eff=eff) if GCl.mockobs.z_snap==z_snap] :
                    List_z = [GCl for GCl in Survey.FilterCluster(minrel=1, eff=feff)] 
                    
                    if len(List_z) > 0:  
                            data_x = [parA(meas).value for  meas in List_z]
                            data_y = [parB(meas).value for  meas in List_z]
                            proxyA,proxyB = parA(List_z[0]), parB(List_z[0]) 
            
                    Data_x.append(data_x)   
                    Data_y.append(data_y)   
                    SurveyHistkwargs.append(Survey.histkwargs)
        
        '''======= New style (begin)'''   
        print(len(Data_x))    
        # Scatter
        for Histkwargs, data_x, data_y in zip(SurveyHistkwargs, Data_x, Data_y):
             
            plotted = axScatter.errorbar( data_x, data_y, fmt = z_symbols[0], alpha=min(3.0*np.power(len(data_x),-0.40),1.))     # , **Survey.scatterkwargs
                            
            try:
                zlist.append(plotted)
            except:
                warnings.warn("%s, efficiency %.2e has no detections." % (Survey.name,eff), DeprecationWarning)
                
            try:
                plotl.append(plotted)
            except:
                warnings.warn("%s, efficiency %.2e has no detections." % (Survey.name,eff), DeprecationWarning)                         
                          

        # Axis
        if proxyA.vrange[0] == proxyA.vrange[1]: proxyA.vrange = [None,None] #Workaround
        if proxyB.vrange[0] == proxyB.vrange[1]: proxyB.vrange = [None,None] #Workaround
        
        if     logplot[0]: proxyA.vrange[0] = 0.1    
        if not logplot[0]: proxyA.vrange[0] = 0.0
                      
        axScatter.set_xlim(left  =proxyA.vrange[0], right=proxyA.vrange[1])
        axScatter.set_ylim(bottom=proxyB.vrange[0],   top=proxyB.vrange[1])  #top=proxyB.vrange[1]  --> but bug
        axScatter.set_ylim(bottom=6e13, top=1e16)   
    
        axScatter.set_xlabel( proxyA.labels(log=log[0]) )
        axScatter.set_ylabel( proxyB.labels(log=log[1]) )
#        plt.tight_layout()
    	    

        lockedaxis = fig.get_axes()[0]
        axScatter.set_autoscale_on(False)
                     
                                                               
        # Histo
        for Histkwargs, data_x, data_y in zip(SurveyHistkwargs, Data_x, Data_y):
                
                
            # Plotted
            print(len(data_x))
            if len(data_x) > 0:
                if logplot[0]:
                    axHistx.hist(data_x, bins=np.logspace( np.log10(axScatter.get_xlim()[0]), np.log10(axScatter.get_xlim()[1]), 10)                     ,  **Histkwargs) 
                    axHistx.hist(data_x, histtype='step', bins=np.logspace( np.log10(axScatter.get_xlim()[0]), np.log10(axScatter.get_xlim()[1]), 10)    ,  color='black', alpha=0.6) 
                    axHistx.set_xscale("log")  #pl.gca().
                else:
                    axHistx.hist(data_x, bins=np.linspace( axScatter.get_xlim()[0], axScatter.get_xlim()[1], 10)                    ,  **Histkwargs) 
                    axHistx.hist(data_x, histtype='step', bins=np.linspace( axScatter.get_xlim()[0], axScatter.get_xlim()[1], 10)  ,  color='black', alpha=0.6)   
                if logplot[1]:   
                    axHisty.hist(data_y, bins=np.logspace( np.log10(axScatter.get_ylim()[0]), np.log10(axScatter.get_ylim()[1]), 10) ,  orientation='horizontal',  weights=[1./max(len(data_y),1)]*len(data_y) ,  **Histkwargs)
                    axHisty.hist(data_y, histtype='step', bins=np.logspace( np.log10(axScatter.get_ylim()[0]), np.log10(axScatter.get_ylim()[1]), 10) ,  orientation='horizontal',  weights=[1./max(len(data_y),1)]*len(data_y) ,  color='black', alpha=0.6)
                    axHisty.set_yscale("log")
                else:
                    axHisty.hist(data_y, bins=np.linspace( axScatter.get_ylim()[0], axScatter.get_ylim()[1], 10) ,  orientation='horizontal',  **Histkwargs)
                    axHisty.hist(data_y, histtype='step', bins=np.linspace( axScatter.get_ylim()[0], axScatter.get_ylim()[1], 10) ,  orientation='horizontal',  color='black', alpha=0.6)
                    
                axHistx.set_xlim(axScatter.get_xlim())
                axHisty.set_ylim(axScatter.get_ylim()) 
        #        
                
                
        '''======= New style (end)'''

#        axScatter.legend(plotl,['CW$_\mathrm{zsnap=0.0}$', 'CW$_\mathrm{zsnap=0.5}$', 'NVSS'], loc=0, frameon=False, handletextpad=0.1,  bbox_to_anchor=(0.25, 0.94), borderaxespad=0.) # loc=2
        axScatter.legend(plotl,['Simu$_\mathrm{all}$', 'Simu$_\mathrm{relicdet}$', 'NVSS'], loc=0, frameon=False, handletextpad=0.1,  bbox_to_anchor=(0.30, 0.74), borderaxespad=0.) # loc=2
    #    red_patch = mpatches.Patch(color='red', label='The red data')
    
    
    
        axScatter.set_autoscale_on(False)
        shape_line_x =  np.linspace(0,1,30)   # [1,2.5,30]
        slope        = 2
        shape_line_y1  =  [  10**(13.8)*(slope**shape_line_x)   for x in shape_line_x]
        axScatter.plot(shape_line_x, shape_line_y1, color='black', alpha=0.7)
        
        # no labels - second try
        axHistx.xaxis.set_major_formatter(nullfmt) # does just remove one label, weird ...
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        #== Save file
        nowfile = 'CraftyPlot_mass_VS_redshift_e%.0f' % (math.log10(eff)*100)     
        #                fig2,ax1,proxyB, lockedaxis = create_Samples_A_histo( SurveySamples, parB, (zrange,colors), eff=eff, **addargs)  
        #                fig.add_subplot(1,1,1)
        #                fig.add_subplot(1,1,1)
        #        plt.figure(fig)
        # ask for the surveys model
        # nowfolder  = Survey.outfolder + '/PostProcessing/%s_%s' % ( parA(meas).dic, parB(meas).dic )   #parA(meas).name, parB(meas).name)
        nowfolder  = SurveySamples[0].outfolder + '/PostProcessing/'  
        iom.check_mkdir(nowfolder)
        print('Gonna save:  %s' % (nowfolder + nowfile))                           
        fig.savefig('%s%s_scatter.png' % (nowfolder,nowfile)) #filename
        fig.savefig('%s%s_scatter.pdf' % (nowfolder,nowfile)) #filename	    #fig.clf()
        
        fig.clf() 

    return 0


def create_Test( SurveySamples, symbolsism, R200exp=False, markers = np.asarray(['.','s']), effi=[], log=[False,False], logplot=[True,True], lockedaxis=False, minrel=1):
    """
    Takes as inputs;
    Surveysample - a list of galaxyclusters
    **args       - optional parameter
    Returns a scatter plot (a fig object) for further modification: linear or logarithmic, in the variables (you have to specify) and saves it under the variable specific name
    
    inspired by : http://www.astrobetter.com/blog/2014/02/10/visualization-fun-with-python-2d-histogram-with-1d-histograms-on-axes/  
             and: http://matplotlib.org/examples/pylab_examples/scatter_hist.html
    """
    
    (zrange,colors,z_symbols) = symbolism
    
    # (r.iner_rat/(r.LAS/(r.dinfo.beam[0]/60.)))
    if len(effi)==0: effi = SurveySamples[0].Rmodel.effList
                 


#    parA = lambda x: x.Mach
#    parA = lambda x: x.T
    plotpairs = [(lambda x: x.Rho   ,  lambda x: x.P_rest, 'Rho_Power'),
                 (lambda x: x.Mach  ,  lambda x: x.P_rest, 'Mach_Power'),
                 (lambda x: x.T     ,  lambda x: x.P_rest, 'T_Power'),
                 (lambda x: x.Area  ,  lambda x: x.P_rest, 'Area_Power'),       
                 (lambda x: x.Area_av,  lambda x: x.P_rest, 'Area_av_Power')
                 ]
    
    cm = plt.cm.get_cmap('RdYlBu')
    
    
    for parA,parB,craftname in plotpairs:

        titlepost = ''


    
        for eff in effi: 
            
            #======= Old style (begin)
    #        fig = plt.figure(facecolor="white", figsize=(5.6, 4.2)) # aspect='equal'
    #        axScatter = fig.add_subplot(1,1,1)
    #        axScatter.tick_params(which='both', direction='in')
    
    
            '''======= New style (begin)'''
            
            nullfmt = NullFormatter()         # no labels
            
            # definitions for the axes
            left, width = 0.1, 0.65
            bottom, height = 0.1, 0.65
            bottom_h = left_h = left + width + 0.02
            
            rect_scatter = [left, bottom, width, height]
            rect_histx   = [left, bottom_h, width, 0.2]
            rect_histy   = [left_h, bottom, 0.2, height]
            
            # start with a rectangular Figure
            fig = plt.figure(1, figsize=(6.5, 6.5))
            
            axScatter = plt.axes(rect_scatter)
            axHistx = plt.axes(rect_histx)
            axHisty = plt.axes(rect_histy)
            
            # ticks inside
            axScatter.tick_params(which='both', direction='in')
            axHistx.tick_params(which='both', direction='in')
            axHisty.tick_params(which='both', direction='in')
            
            
            #show grid
            axScatter.grid()
            axHistx.grid()
            axHisty.grid()
            
            # no labels
            axHistx.xaxis.set_major_formatter(nullfmt)
            axHisty.yaxis.set_major_formatter(nullfmt)
            
            '''======= New style (end)'''
    #        axScatter.set_xlim( (0,0.3) )
    #        ax.set_ylim( (2e-2,1) )
    #        ax.set_autoscale_on(False)
        
            if logplot[0]: axScatter.set_xscale('log')
            if logplot[1]: axScatter.set_yscale('log')
        
    
            plotl  = []
            zlist  = []
            Data_x = []
            Data_y = []
            SurveyHistkwargs = []
        #   sc = plt.scatter(xy, xy, c=z, vmin=0, vmax=20, s=35, cmap=cm)
        #   cbar = plt.colorbar(sc)
        #   cbar.solids.set(alpha=1)
        
            for Survey in SurveySamples:
                
                print('%s::Survey.name:' % (craftname,Survey.name))
                
                if Survey.Rmodel.simu: 
                    feff = eff
                    scale = min(effi)/eff   
                else:
                    feff = 1.
                    scale = 1.
                if Survey.name == 'NVSS':
                    Survey.scatterkwargs.update( color='black', ecolor='black', alpha=0.4, fmt=None)
                    Survey.histkwargs.update( color='black', alpha=0.4)
                    del Survey.scatterkwargs['color']
                    del Survey.histkwargs['color']
                else:
                    Survey.scatterkwargs.update( color='r', ecolor='r', fmt=None)
                    Survey.histkwargs.update( color='r')
                    del Survey.scatterkwargs['color']
                    del Survey.histkwargs['color']
                del Survey.scatterkwargs['fmt']
                    # make the right colors
    
                   
                data_x = []; data_y=[]
    #            List_z = [GCl for GCl.filterRelics(eff=eff) in l for GCl in Survey.FilterCluster(minrel=minrel, eff=feff)] 
                List_zz = []
                for GCl in Survey.GCls:
                        List_zz.append([relic for relic in GCl.filterRelics(eff=feff) if (relic.region.rtype.classi != 0)]  )
                List_z = [item for sublist in List_zz for item in sublist]
                print(Survey.name, feff, len(List_z))
                
                if len(List_z) > 0:  
                        data_x = [parA(meas).value        for  meas in List_z]
                        data_y = [parB(meas).value*scale  for  meas in List_z]
                        proxyA,proxyB = parA(List_z[0]), parB(List_z[0]) 
        
                Data_x.append(data_x)   
                Data_y.append(data_y)   
                SurveyHistkwargs.append(Survey.histkwargs)
                  
            
            '''======= New style (begin)'''   
                                 
                             
            for Histkwargs, data_x, data_y in zip(SurveyHistkwargs, Data_x, Data_y):
            # Scatter
                
                if len(data_x) > 0:
                    plotted = axScatter.errorbar( data_x, data_y, fmt = z_symbols[0], alpha=min(3.0*np.power(len(data_x),-0.40),1.))     #, **Survey.scatterkwargs
                                    
                    try:
                        zlist.append(plotted)
                    except:
                        warnings.warn("%s, efficiency %.2e has no detections." % (Survey.name,eff), DeprecationWarning)
                        
                    try:
                        plotl.append(plotted)
                    except:
                        warnings.warn("%s, efficiency %.2e has no detections." % (Survey.name,eff), DeprecationWarning)           
                        
                        
            #    red_patch = mpatches.Patch(color='red', label='The red data')
        
            # This we have to fix, for many measurands [1,1] is given instead of [None,None]
            if proxyA.vrange[0] == proxyA.vrange[1]: proxyA.vrange = [None,None] #Workaround
            if proxyB.vrange[0] == proxyB.vrange[1]: proxyB.vrange = [None,None] #Workaround
            
                       
            axScatter.set_xlim(left  =proxyA.vrange[0], right=proxyA.vrange[1])
            axScatter.set_ylim(bottom=proxyB.vrange[0],   top=proxyB.vrange[1])  #top=proxyB.vrange[1]  --> but bug
    
        
            axScatter.set_xlabel( proxyA.labels(log=log[0]) )
            axScatter.set_ylabel( proxyB.labels(log=log[1]) )
    #        plt.tight_layout()
        	    
            lockedaxis = fig.get_axes()[0]
            
            
            axScatter.set_autoscale_on(False)    
            
            '''  Histo'''                                     
            for Histkwargs, data_x, data_y in zip(SurveyHistkwargs, Data_x, Data_y):
    
                
                if len(data_x) > 0:                                            
    
                    axHistx.hist(data_x, bins=np.logspace( np.log10(axScatter.get_xlim()[0]), np.log10(axScatter.get_xlim()[1]), 10)                            ,  weights=[1./max(len(data_x),1)]*len(data_x) ,  **Histkwargs) 
                    axHisty.hist(data_y, bins=np.logspace( np.log10(axScatter.get_ylim()[0]), np.log10(axScatter.get_ylim()[1]), 10) ,  orientation='horizontal',  weights=[1./max(len(data_y),1)]*len(data_y) ,  **Histkwargs)
                    
                    axHistx.hist(data_x, histtype='step', bins=np.logspace( np.log10(axScatter.get_xlim()[0]), np.log10(axScatter.get_xlim()[1]), 10)                            ,  weights=[1./max(len(data_x),1)]*len(data_x) ,  color='black', alpha=0.6) 
                    axHisty.hist(data_y, histtype='step', bins=np.logspace( np.log10(axScatter.get_ylim()[0]), np.log10(axScatter.get_ylim()[1]), 10) ,  orientation='horizontal',  weights=[1./max(len(data_y),1)]*len(data_y) ,  color='black', alpha=0.6)
        
                    axHistx.set_xlim(axScatter.get_xlim())
                    axHisty.set_ylim(axScatter.get_ylim())
            #        
                    axHistx.set_xscale("log")  #pl.gca().
                    axHisty.set_yscale("log")
            '''======= New style (end)'''
    
            fig.legend(plotl,['CW', 'NVSS'], loc=0, frameon=False, handletextpad=0.1,  bbox_to_anchor=(0.24, 0.73), borderaxespad=0.) # loc=2
            
       
            # no labels - second try
            axHistx.xaxis.set_major_formatter(nullfmt) # does just remove one label, weird ...
            axHisty.yaxis.set_major_formatter(nullfmt)
            
            #== Save file
            nowfile = 'CraftyPlot_%s_e%.0f' % (craftname, math.log10(eff)*100)     
            #                fig2,ax1,proxyB, lockedaxis = create_Samples_A_histo( SurveySamples, parB, (zrange,colors), eff=eff, **addargs)  
            #                fig.add_subplot(1,1,1)
            #                fig.add_subplot(1,1,1)
            #        plt.figure(fig)
            # ask for the surveys model
            # nowfolder  = Survey.outfolder + '/PostProcessing/%s_%s' % ( parA(meas).dic, parB(meas).dic )   #parA(meas).name, parB(meas).name)
            nowfolder  = SurveySamples[0].outfolder + '/PostProcessing/'  
            iom.check_mkdir(nowfolder)
            print('Gonna save:  %s' % (nowfolder + nowfile) )                          
            fig.savefig('%s%s_scatter.png' % (nowfolder,nowfile)) #filename
            fig.savefig('%s%s_scatter.pdf' % (nowfolder,nowfile)) #filename	    #fig.clf()
            
            fig.clf() 
            
  
'''=== Section of single object analysis ==='''   
'''========================================='''
    

def joinpandas(pdframes):
    pdframe_combined = None
    
    for pdframe in pdframes:
        if pdframe_combined is None:
            pdframe_combined  = pdframe
        else:
            pdframe_combined = pdframe_combined.append(pdframe)

    return pdframe_combined
    
def create_scattermatrix( SurveySamples, plotmeasures, suffix=''):
     
    ''' Creates a scatter matrix, off a list of quantities ... nice! 
    Input: SurveySamples ... there is currently no differnciation between different Survey Samples (symbolwise or else)
    ''' 
    from pandas.tools.plotting import scatter_matrix
 
    pdframes = [survey.fetchpandas(plotmeasures) for survey in SurveySamples]
    pdframe_combined = joinpandas(pdframes)                    

    print(len(SurveySamples))    
           
            
    ''' Examples of additional plots 
    def hexbin(x, y, color, **kwargs):
         cmap = sns.light_palette(color, as_cmap=True)
         plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)
        
    def corplot(x, y, color, **kwargs):
         cmap = sns.light_palette(color, as_cmap=True)
         plt.imshow(np.abs([x,y].corr()), cmap=cmap, **kwargs) 
    '''
    
    

    import seaborn as sns; sns.set(style="ticks", color_codes=True)
    from itertools import cycle
    
    pdframe_combined.to_csv(path_or_buf='/data/Test-%s.csv' % (SurveySamples[0].name))
    print(pdframe_combined.Survey.unique())
    g = sns.PairGrid(pdframe_combined, hue="Survey", palette="Set2",dropna=True)     
    g = g.map_upper (sns.regplot, scatter_kws={'edgecolors':"white","linewidth":1,"alpha":0.3})  #plt.scatter , , edgecolor="white"
    g = g.map_diag(sns.kdeplot, lw=3, legend=False, alpha=0.5)  #histtype="step"  {'cmap':['Blues_d','Blues']}
    
    # Taken from https://stackoverflow.com/questions/40726733/plotting-multiple-datasets-on-a-seaborn-pairgrid-as-kdeplots-with-different-colo
    def make_kde(*args, **kwargs):    
        sns.kdeplot(*args, cmap=next(make_kde.cmap_cycle), **kwargs)
        
    pdframe_combined.Survey.unique()    
    colorsmaps = ('BuGn','Oranges', 'Red') 
        
    make_kde.cmap_cycle = cycle(colorsmaps[0:len(pdframe_combined.Survey.unique())])    #, 'Reds_r'
        
    g = g.map_lower(make_kde, alpha=0.5) #cmap="Blues", shade=True,   color=...
   

    
    #== Save file
    nowfile = 'Scattermatrix'
    #                fig2,ax1,proxyB, lockedaxis = create_Samples_A_histo( SurveySamples, parB, (zrange,colors), eff=eff, **addargs)  
    #                fig.add_subplot(1,1,1)
    #                fig.add_subplot(1,1,1)
    #        plt.figure(fig)
    # ask for the surveys model
    # nowfolder  = Survey.outfolder + '/PostProcessing/%s_%s' % ( parA(meas).dic, parB(meas).dic )   #parA(meas).name, parB(meas).name)
    nowfolder  = SurveySamples[-1].outfolder + '/PostProcessing/'  
    iom.check_mkdir(nowfolder)
    print('Gonna save:  %s' % (nowfolder + nowfile)  )                         
    plt.savefig('%s%s%s.png' % (nowfolder,nowfile,suffix),dpi=400) #filename
    plt.savefig('%s%s%s.pdf' % (nowfolder,nowfile,suffix)) #filename	    #fig.clf()

    plt.clf() 
#    fig.clf() 
            

#
#def oldscattermatrix(R200exp=False, markers = np.asarray(['.','s']), log=[False,False], logplot=[True,True], lockedaxis=False, minrel=1):
#    '''    
#        if 1 == 2:
#            
#            if cc == 0:    
#                fig, axes = scatter_matrix(pdframe, alpha=0.65, figsize=(12, 12), c=col, cmap=massmap, edgecolors=(0,0,0,0.3)) #    ps=6,  , diagonal='kde', , **{'scale':'log'}
#        
#                if masscolor:
#                    cbarplot = plt.scatter(cvalues_mass, cvalues_mass, c=cvalues_mass, cmap=massmap, vmin=0.0, vmax=1.0) #viridis
#                    colorbar_ax = fig.add_axes([0.11, 0.03, 0.80, 0.03])  # [left, bottom, width, height]
#                    cb =  fig.colorbar(cbarplot, cax=colorbar_ax, orientation="horizontal") 
#                    cb.set_label('$\log_{10} (M_\mathrm{cl}/M_\odot)$' )
#            
#                a = np.linspace(0.0,1.0,0.1, endpoint=True)
#                b = np.ones_like(a)
#                cbarplot2 = plt.scatter(a, a,  c=b, cmap=cormap, vmin=0.0, vmax=1.0)
#                colorbar_ax2 = fig.add_axes([0.11, 0.93, 0.80, 0.03])  # [left, bottom, width, height]
#                cb =  fig.colorbar(cbarplot2, cax=colorbar_ax2, orientation="horizontal")
#                cb.set_label('Correlation strength')
#            else:
#                scatter_matrix(pdframe, alpha=0.65, figsize=(12, 12), c=col, cmap=massmap, edgecolors=(0,0,0,0.3))
#
#            plt.subplot_tool()
#            plt.show()
#    '''
#    cm = plt.cm #.get_cmap('RdYlBu')
##     titlepost = ''
#    
#    colList = ['b','r','g','greyy','yellow']
#
#    '''======= New style (begin)'''   
#    masscolor = False
#    
#    if masscolor:
#        col =  cvalues_mass #cm.plasma    
#    else: 
#        col       = colList[cc]
#    
#    cormap    = cm.magma #cm.viridis
#    massmap   = cm.plasma #cm.magma    #cm.viridis
#
#     
#    ''' A workaround, you optionally could initialize the scatterplot with the log argument '''
#    for i, axs in enumerate(axes):
#        for j, ax in enumerate(axs):
#            if i < j and cc==0:  # only the scatter plots
##                    ax.set_xscale('log')
##                    ax.set_yscale('log')
#
##                    fig.delaxes(ax)
#
#                ax.clear()
#                ''' Plot correlation '''
#                ax.imshow([[np.abs(pdframe.corr().iloc[i,j])]], vmin=0.0, vmax=1.0, cmap=cormap)
#                ax.text(0.5, 0.5,'%.2f' % (np.abs(pdframe.corr().iloc[i,j])), horizontalalignment='center',verticalalignment='center',transform=ax.transAxes, fontsize=30, color='white')
#
#            if i == j and cc>0:  # only the scatter plots
#                ax.hist(pdframe[pdframe.index == i], normed=1, facecolor=colList[cc], alpha=0.6)
#                
#            if i > j and cc >0:       
#                ax.scatter(pdframe[pdframe.index == i], pdframe[pdframe.index == j], color=colList[cc], alpha=0.6)

