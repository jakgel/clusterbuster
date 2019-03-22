#!/usr/bin/env python

""" 
Created on 2017
@author: jakobg

This .py provides the functionalities of higher-level data products from the pickled relic catalogues, including:
- Creating pandas data tables (implement saving them as .csv or o.ods files)
- Creating .pandas scatter matrixes
- Creating .fits or .png images of the simulated objects and their galaxy clusters
- Creating ...
"""

from __future__ import division,print_function

import copy
import os
import warnings
import aplpy
import collections  # to assert if we deal with an list of surveys or a single object

import clusterbuster.surveyclasses as cbclass
import clusterbuster.dbclasses as dbc
import clusterbuster.iout.misc as iom
import clusterbuster.maput as maput
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np


from matplotlib import cm
from matplotlib.pyplot import tight_layout
from matplotlib.ticker import NullFormatter
from matplotlib import rc
from astropy import units as u
import matplotlib.patches as patches
from matplotlib import transforms as mtransforms

import surveysim.music2.mockobsxray as Xray
from PyPDF2 import PdfFileMerger
from scipy import stats

import seaborn as sns;
from itertools import cycle
from pandas.tools.plotting import scatter_matrix
from scipy.stats import norm

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# ============== Reads relic information out of an ds9 region file
def readDS9relics(regfile, spixel, center, pixref, Test=False):
    contours, contourinfos = iom.readDS9regions(regfile, spixel, center, pixref)
    contoursWCS, _ = iom.readDS9regions(regfile, spixel, center, pixref, pixelcoords=False)

    rinfo = []
    for ii, info in enumerate(contourinfos):
        # info       =  info.split(', ')

        try:
            info = info.split(', ')  # get infolist
            info[2] = info[2].replace('alpha ', '')  # remove this alpha string
            split = info[2].split(' ')  # split alpha into substrings
            if len(split) > 1:
                alpha = split[0]
                alpha_err = split[1]
            else:
                alpha = split[0]
                alpha_err = 0

            reg = cbclass.RelicRegion(name=info[0], cnt=[contours[ii]], cnt_WCS=[contoursWCS[ii]], rtype=int(info[1]),
                                      alpha=alpha, alphaFLAG=('false' not in alpha.lower()), alpha_err=alpha_err)
        except:
            reg = cbclass.RelicRegion(name='', cnt=[contours[ii]], cnt_WCS=[contoursWCS[ii]], rtype=-1, alphaFLAG=False)

        if ('test' not in info[0].lower()) or not Test:
            rinfo.append(reg)

    return rinfo


def plot_RelicEmission_polar(surveys, compsurvey=None, single=False, modeltext=True, additive=False,
                             aligned=False, cbar=True, addinfo=False, mirrored=False, plottype='flux',
                             title="Polar binned radio relic flux", dpi=None, add_pi=1/2,
                             Histo=dbc.Histogram2D(), suffix='', conv=0.127):
    """ possible inprovements: http://stackoverflow.com/questions/22562364/circular-histogram-for-python
    minrel : minimal number of relics to be consided for the histogram
    addpi: Additional Rotation;  added anticlockwise! We want to turn by 90 degree anticlockwise

    if surveys is a list of surveys this function will plot averaged quantities
    """
    plt.style.use('default')
    if not isinstance(surveys, collections.Iterable):
        surveys = [surveys]

    if surveys[0].hist_main is not None:
        Histo = surveys[0].hist_main

    """ Plotting and normalization of combined histogram"""
    expPlot  = 0.45
    dist_text_params = 0.08
    cmap     = cm.viridis
    yticks   = [0.4, 0.8, 1.2]
    addangle = int(aligned)*np.pi*add_pi

    halfHists = []
    radials   = []
    stats     = []
    outformats = ['pdf', 'png']


    if compsurvey is not None:
        _ , (comprad, _), comppol, _, _ = compsurvey.polar(conv=conv)
        deviations = []

    for survey in surveys:

        halfes={'First': Histo.bins[0][int(len(Histo.bins[0]) / 2):],
                'Second': Histo.bins[0][0:int(len(Histo.bins[0]) / 2) + 1]}
        if mirrored:
            half_main = halfes['First']
            half_seco = halfes['Second']
        else:
            half_main = halfes['Second']
            half_seco = halfes['First']

        survey.set_binning(Histo)
        nowfolder = '%s/Relics_polar/' % survey.outfolder
        iom.check_mkdir(nowfolder)

        buckedfolder = os.path.abspath(os.path.join(survey.outfolder, '..', 'bucket'))
        iom.check_mkdir(buckedfolder)
        if single:
            for ii, GCl in enumerate(survey.FilterCluster()):

                """ Deriving the histogram should be a functionality of the survey or the relic cluster, so this should become outdated 
                Beware that at this point, survey.Hist and gcl.Hist are the same objects!
                
                """
                GCl.updateInformation(Filter=True)
                if GCl.histo is not None and np.sum(GCl.histo.hist) != 0:
                    inner = Histo.bins[1][0:-1]
                    outer = Histo.bins[1][1::]
                    angle = Histo.ticks[0]
                    angles, z_inner = np.meshgrid(angle, inner, sparse=True)
                    angles, z_outer = np.meshgrid(angle, outer, sparse=True)
                    shiftHist = np.roll(GCl.histo.hist.T, -int(aligned*(GCl.relic_pro_index)), axis=1) / survey.AreaHist**(survey.expA)  #+1e-10
                    # Plots the single clusters
                    fig, ax = plt.subplots(figsize=(14,14), subplot_kw=dict(projection='polar'), dpi=dpi)
                    ax.pcolormesh(Histo.bins[0], Histo.bins[1],  shiftHist, cmap=cmap)
                    ax.set_theta_offset(addangle)
                    ax.annotate("", xy=(int(not aligned)*(GCl.relic_pro_angle), 1.5), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))
                    ax.arrow(0, 0, 0, 0.5, linewidth=3, width=0.005, transform=mtransforms.Affine2D().translate(int(not aligned)*(GCl.relic_pro_angle), 0) + ax.transData)
                    ax.text(0.01, 1.05, '%s' % (GCl.name.replace('_', ' ')), fontsize=20, transform=ax.transAxes)
                    if addinfo:
                        ax.text(0.3, 0.9, 'Summed relic flux: %.2e Jy' % (np.sum(Histo.hist)), fontsize=20, transform=ax.transAxes, color='w')
                        ax.text(0.3, 0.87, 'Ratio pro: %.2e  anti: %.2e' % (GCl.ratio_pro(), GCl.ratio_anti()), fontsize=20, transform=ax.transAxes, color='w')
                    if title is not None:
                        ax.set_title(title, va='bottom')
                    ax.set_rticks(yticks)
                    ax.tick_params(axis='x', labelsize=25)
                    ax.tick_params(axis='y', colors='white', labelsize=25, pad=23)
                    ax.set_rlabel_position(89.9)
                    ax.grid(True)
                    tight_layout()
                    for ftype in outformats:
                        plt.savefig('%s%s-polar%s.%s' % (nowfolder, GCl.name, suffix, ftype))

                    fig.clf()

        if additive:
            _, (radial, radial_tickz), halfHist, stat, mesh = survey.polar(aligned=True, mirrored=mirrored, mode=plottype, conv=conv)

            if halfHist is not None:
                fig, ax = plt.subplots(figsize=(14, 14), subplot_kw=dict(projection='polar'), dpi=dpi)  #,subplot_kw=dict(projection='polar')
                near_max = np.percentile(halfHist, 99)
                if mirrored:
                    meshed = ax.pcolormesh(Histo.bins[0][int(len(Histo.bins[0])/2):], Histo.bins[1], halfHist, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=near_max)
                    if compsurvey:
                        ax.pcolormesh(Histo.bins[0][0:int(len(Histo.bins[0])/2)+1], Histo.bins[1], comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=near_max)
                else:
                    meshed = ax.pcolormesh(Histo.bins[0][0:int(len(Histo.bins[0])/2)+1], Histo.bins[1], halfHist, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=near_max) #, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-1.0, vmax=1.0),cmap='RdBu_r'
                    if compsurvey:
                        ax.pcolormesh(Histo.bins[0][int(len(Histo.bins[0])/2):], Histo.bins[1], comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=near_max)


                #meshed = ax.pcolormesh(half_main, Histo.bins[1], halfHist, cmap=cmap,
                #                       norm=colors.PowerNorm(gamma=expPlot), vmax=near_max)
                #if compsurvey is None:
                #    ax.set_thetamin(0)
                #    ax.set_thetamax(180)
                #else:
                #    ax.pcolormesh(half_seco, Histo.bins[1], comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot),
                #                  vmax=near_max)

                ax.set_theta_offset(addangle)
                ax.set_rticks(yticks)
                ax.tick_params(axis='x',                 labelsize=25)
                ax.tick_params(axis='y', colors='white', labelsize=25, pad=23)
                ax.set_rlabel_position(89.9)
                ax.grid(True)
                tight_layout()
                if addinfo:
                    ax.text(0.3, 0.9, 'Summed relic flux: %.2e Jy (all cluster)' % stat, fontsize=20, transform=ax.transAxes, color='w')
                if cbar:
                    fig.colorbar(meshed)
                if title is not None:
                    ax.set_title(title, va='bottom')
                for ftype in ['pdf', 'jpg']:
                    plt.savefig('%s/%s-polar%s.%s' % (nowfolder, survey.name, suffix, ftype))
                fig.clf()

            halfHists.append(halfHist)
            radials.append(radial)
            stats.append(stat)
            if compsurvey is not None:
                deviations.append(np.sum(np.abs(radial-comprad)))

            # plot ratio of relics flux,
            # plot average/median pro relic distance
            # plot sigma pro rleic distance
            # plot skew

    """ Colorsheme from seaborn """
#    cmap = ListedColormap(sns.color_palette('deep'))
    if len(radials) > 0:

        """ Radial plot """
        scale = 0.8
        plt.rcParams['figure.figsize'] = [7 * scale, 4 * scale]
        plt.subplots_adjust(hspace=0.4)


        scaling=1
        fig, (ax1) = plt.subplots(1, 1, figsize=(7*scaling,4.2*scaling), dpi=dpi)
        ax1.set_xlabel('Distance [$R_{200}$] along pro-relic axis')
        ax1.set_ylabel('Weighted signal S')

        if compsurvey is not None:
            ax1.plot(radial_tickz, comprad, alpha=0.6, color='blue')
            ax1.fill_between(radial_tickz, comprad, color="blue", alpha=0.2)

            patch1 = patches.Patch(color='blue', alpha=0.2, label=compsurvey.name_short)
            patch2 = patches.Patch(color='red', alpha=0.2, label=survey.name_short)
            plt.legend(handles=[patch1,patch2])

        for radial in radials:
            ax1.plot(radial_tickz, radial, alpha=np.sqrt(1/len(radials)), c='grey')   # color=cmap(0.0)
            ax1.fill_between(radial_tickz, radial, color="red", alpha=0.2)

        #ax1.set_xticks(yticks + [0] + [-y for y in yticks])
        ax1.set_xlim(-Histo.bins[1][-1],Histo.bins[1][-1])
        ax1.set_ylim(bottom=0)
        ax1.tick_params(direction="in", which='both')
        ax1.tick_params(direction="in", which='major', right=False, top=True, labelright=False)

        """ Textlabels """
        if modeltext and survey.Rmodel.simu:
            mod = survey.Rmodel
            kwargs = {'verticalalignment':'bottom', 'horizontalalignment':'right', 'transform':ax1.transAxes, 'color':'black', 'fontsize':12, 'alpha':0.8}
            ax1.text(0.35, 0.90, '$\log_{10}(eff) =%+.2f$,' % (np.log10(mod.effList[0])),   **kwargs)
            ax1.text(0.35, 0.90-1*dist_text_params, '$\log_{10}(B_0) =%+.2f$,' % (np.log10(mod.B0)),**kwargs)
            ax1.text(0.35, 0.90-2*dist_text_params, '$\kappa       = %+0.2f$,' % (mod.kappa),       **kwargs)
            if isinstance(mod, cbclass.PreModel_Hoeft):
                ax1.text(0.35, 0.90-4*dist_text_params, '$t_{1;2}  = %0.3f\,;\,%0.3f$,'% (mod.t0, mod.t1), **kwargs)
                ax1.text(0.35, 0.90-5*dist_text_params, 'ratio$\\mathrm{_{pre}}  = %.1e$,' % mod.ratio, **kwargs)
            if compsurvey is not None:
                ax1.text(0.35, 0.90-3*dist_text_params, '$\Delta\\,\\mathrm{signal}  = %0.2f$ ' % (np.average(deviations)), **kwargs)

            if survey.Rmodel.pre:
                """ NOT implemented yet """
#                print( p0, pre, p_sigma, sigmoid_0, sigmoid_width )
    #            ax2.set_yscale('log')

        for ftype in outformats:
            plt.savefig('%s/%s-sumprofile%s.%s' % (nowfolder, survey.name, suffix, ftype))
            plt.savefig('%s/%s-sumprofile%s.%s' % (buckedfolder, survey.name, suffix, ftype))


        # Statistics, contribution
        fig, (ax1) = plt.subplots(1, 1, figsize=(7, 4), dpi=dpi)

        for stat in stats:
            statistic = np.divide(stat, np.sum(stat))
            ax1.hist(statistic, color='b', alpha=1/len(stats), bins='auto')  # arguments are passed to np.histogram
        for ftype in outformats:
            plt.savefig('%s/%s-sumsstats%s.%s' % (nowfolder, survey.name, suffix, ftype))
            plt.savefig('%s/%s-sumsstats%s.%s' % (buckedfolder, survey.name, suffix, ftype))


        """ Polar plot"""
        halfHist = np.sum(halfHists, axis=0)/len(halfHists)
        fig, ax = plt.subplots(figsize=(14,14), subplot_kw=dict(projection='polar'), dpi=dpi)  #,subplot_kw=dict(projection='polar')
        near_max = np.percentile(halfHist, 99)
        meshed = ax.pcolormesh(half_main, Histo.bins[1], halfHist, cmap=cmap,
                               norm=colors.PowerNorm(gamma=expPlot), vmax=near_max)
        if compsurvey is None:
            ax.set_thetamin(0)
            ax.set_thetamax(180)
        else:
            ax.pcolormesh(half_seco, Histo.bins[1], comppol, cmap=cmap, norm=colors.PowerNorm(gamma=expPlot), vmax=near_max)

        ax.set_theta_offset(addangle)
        ax.set_rticks(yticks)
        ax.tick_params(axis='x',                 labelsize=25, pad=23)
        ax.tick_params(axis='y', colors='white', labelsize=25, pad=23)
        ax.set_rlabel_position(89.9)
        ax.grid(True)
        tight_layout()
        if addinfo:
            ax.text(0.3, 0.9, 'Summed relic flux: %.2e Jy (all cluster)' % stat, fontsize=20, transform=ax.transAxes, color='w')
        if cbar:
            fig.colorbar(meshed)
        if title is not None:
            ax.set_title(title, va='bottom')
        for ftype in outformats:
            plt.savefig('%s/%s-polar-sums%s.%s' % (nowfolder, survey.name, suffix, ftype))
            plt.savefig('%s/%s-polar-sums%s.%s' % (buckedfolder, survey.name, suffix, ftype))

        fig.clf()


    """From https://stackoverflow.com/questions/51871420/matplotlib-polar-histogram-has-shifted-bins/51876418
    Includes the binning and rotation of pixelized images ... could be doable with the to fits function"""
    if 1 == 2:
        plt.clf()
        bins_number = 10
        width = 2 * np.pi / bins_number
        ax = plt.subplot(1, 1, 1, projection='polar')
        bars = ax.bar(bins[:bins_number], n, width=width, bottom=0.0)
        for bar in bars:
            bar.set_alpha(0.5)
        plt.show()

    return 0

def plot_Clusters_subfigure(survey, ids):
    pass
    """
    fig = plt.figure(figsize=(8, 10))
    gc = []

    import copy
    survey_use = copy.deepcopy(survey)

    #for GCl in survey_use.GCls:
    #    if GCl.name == ids:
    #        gc.append(aplpy.FITSFigure(GCl.mapdic['Diffuse'], subplot=[0.05, 0.05, 0.9, 0.3], figure=fig))

    path = "/data/ClusterBuster-Output/NVSS/Images"

    gc.append(aplpy.FITSFigure(img1, subplot=[0.05, 0.05, 0.9, 0.3], figure=fig))
    gc.append(aplpy.FITSFigure(img2, subplot=[0.05, 0.35, 0.9, 0.3], figure=fig))

    for subfigure in gc:
        subfigure.recenter(ra, dec, radius=0.5)
        subfigure.tick_labels.hide()
        subfigure.axis_labels.hide()
    """


def plot_Clusters(survey, dynamicscale=False, subtracted=True, relicregions=False, DS9regions=False, diamF=2.6,
                  colorbar=False, beam=True, shapes=False, recenter=True, infolabel = False, sectors=False,
                  xray=False, highres=False, show_rot=False, vectors=False, label_sheme='balanced',
                  filterargs={'zborder': 0, 'ztype': '>', 'minimumLAS': 4, 'GClflux': 20, 'index': None}):
    print('plot_Clusters:BUGS: THere not always the full contours shown')

    sns.set(style="white", color_codes=True)
    pdfs = []
    laargs = {'color': '#DDDDDD'}      # line arguments
    ciargs = {'color': '#DDDDDD'}      # arguments for the circle/centre area
    baargs = {'color': '#DDDDDD'}      # argument for the scale bar
    cmap = 'afmhot'
    cnt_color = 'green'
    cnt_color_sub = 'red'
    color_offset = [0,0]
    if label_sheme == 'dark':
        laargs.update({'color': 'black'})
        ciargs.update({'color': '#111111'})
        baargs.update({'color': '#111111'})
    elif label_sheme == 'bright':
        laargs.update({'color': 'w'})
        ciargs.update({'color': 'w'})
        baargs.update({'color': 'w'})
    elif label_sheme == 'balanced':
        laargs.update({'color': 'snow'})
        ciargs.update({'color': 'snow'})
        baargs.update({'color': 'snow'})
        cmap = "cubehelix" #"gnuplot"
        cnt_color = "honeydew" #greenyellow
        cnt_color_sub = 'darkred' #'salmon'
        color_offset = [1,4]


    for GCl in survey.FilterCluster(**filterargs):

        # Outdated? Because else changes will influence all galaxy clusters, you knwo class referencing in python, i.e.
        GCl = copy.deepcopy(GCl)
        if xray:
            """X-Ray"""
            if 'Brems' not in GCl.mapdic:
                savefolder = survey.outfolder
                Xray.Run_MockObs_XRay(GCl, savefolder, verbose=False)
                """ Is the mapdic updated ? """
                GCl.mapdic['Background'] = GCl.mapdic['Brems']
        else:
            GCl.mapdic['Background'] = GCl.mapdic['Diffuse']

        #=== Sets some plotting parameters
        kpc = GCl.cosmoPS*3600  # kiloparsec per degree

        R200 = GCl.R200()
        if R200 <= 0:
            print('For GCl %s no mass proxy and hence no virial radius is known. For plottting issues we set R200=1600 kpc.' % (GCl.name))
            R200 = 1600
        radius = R200/kpc
        diam = diamF*R200/kpc

        f = aplpy.FITSFigure(GCl.mapdic['Background']) #dimensions=[0, 1],, slices=[10, 10], , slices=[10, 10]

        if recenter:
            print('Recentering', GCl.name, GCl.RA(),GCl.Dec(),diam)
            f.recenter(GCl.RA(), GCl.Dec(), width=diam, height=diam) # radius is also possible!

        if survey.Rmodel.simu:
            f.axis_labels.set_xtext('Coordinate 1')
            f.axis_labels.set_ytext('Coordinate 2')
            f.tick_labels.hide()
        else:
            f.tick_labels.hide()
            #f.tick_labels.set_xformat("dd:mm:ss")
            #f.tick_labels.set_yformat("dd:mm:ss")
            f.axis_labels.hide()

        # The basic image
        if dynamicscale:
            vmax = np.max(f._data)
            levels = [GCl.dinfo.limit*l for l in [survey.m_cnt**n for n in np.arange(0,16)]]  #0,16
        else:
            vmax = survey.emi_max
            levels = survey.cnt_levels

        vmin = 0.6 * GCl.dinfo.rms   #0.25
        vmid = -2   #-2   #
        exponent = np.log(vmax/vmin)

        if highres:
            key = "Raw"
            key_comp = "CompModell"
            levels = [levels[0] / 8, levels[0] / 4, levels[0] / 2] + levels
        else:
            key = "Diffuse"
            key_comp = "Subtracted"



        if key_comp in GCl.mapdic and subtracted:
            f.show_contour(GCl.mapdic[key_comp], linewidth=0.15, overlap=True, levels=levels, colors=cnt_color_sub)


        if not xray:
            cbar_text = 'flux density in [Jy/beam}'

            for relic in GCl.filterRelics():
                pixelcnt = np.transpose(np.squeeze(relic.cnt))
                wcscnts = f.pixel2world(pixelcnt[0,:], pixelcnt[1,:])
                wcscnts = np.asarray([ (x,y) for x,y in zip(wcscnts[0],wcscnts[0]) ]).T
                f.show_polygons([wcscnts], lw=2, color='white')
            addargs = {'vmid': vmid, 'vmin': vmin, 'vmax': vmax, 'stretch': 'log', 'exponent': exponent}

            """ It seems like you can only have one interactive contours """
            print(vmin,vmid,vmax)
            f.show_colorscale(vmid=vmid, vmin=vmin, vmax=vmax,  stretch='log', exponent=exponent, cmap=cmap)
            print(levels, survey.cnt_levels)
            f.show_contour(GCl.mapdic['Diffuse'], linewidth=0.15, overlap=True, levels=levels, colors=cnt_color, smooth=1)
        else:
            cbar_text = '$\log_{10}(P_\\mathrm{Brems,bol}$ in restframe) [arbitrary unit]'
            if 'MUSIC' in survey.name or 'Threehundret' or 'ShockTest' in survey.name:
                vmin_xr = 2.5+color_offset[0]
                vmax_xr = 9.7+color_offset[1]     #6.2
                vmid_xr = -1.5
            else:
                vmin_xr = -2+color_offset[0]
                vmax_xr = 5.+color_offset[1]
                vmid_xr = -4.5

            exponent = np.log(max(vmax/vmin, 1.0001))

            f.show_colorscale(vmid=vmid_xr, vmin=vmin_xr, vmax=vmax_xr,  stretch='log', exponent=exponent, cmap=cmap)  #gist_heat
            if key in GCl.mapdic:
                f.show_contour(GCl.mapdic[key], linewidth=0.15, overlap=True, levels=levels, colors=cnt_color)

        if shapes:
            for relic in GCl.filterRelics(maxcomp=100):
                vlen = (np.sqrt(relic.iner_rat())*relic.LLS + 0.05*R200)/kpc
                f.show_arrows(relic.RA(), relic.Dec(), -relic.eigvecs[1][0]*vlen, relic.eigvecs[1][1]*vlen) #-np.cos(relic.Dec*np.pi/180)*
                f.add_label(relic.RA-relic.eigvecs[1][0]*vlen, relic.Dec+relic.eigvecs[1][1]*vlen, '$s = %.2f$' % (relic.iner_rat()), size='x-large', **laargs)

        # The Jakobs circle OR (virial) radius
        f.show_circles(GCl.RA(), GCl.Dec(), radius, linestyle='--', **ciargs)
        if sectors:
            GCl.relics_polarDistribution(histo=survey.hist_main)
            P1 = [GCl.RA() - np.cos(GCl.relic_pro_angle)*radius*1.05/np.cos(np.radians(GCl.Dec())), GCl.Dec() + np.sin(GCl.relic_pro_angle)*radius*1.0]
            P2 = [GCl.RA() + np.cos(GCl.relic_pro_angle)*radius*1.05/np.cos(np.radians(GCl.Dec())), GCl.Dec() - np.sin(GCl.relic_pro_angle)*radius*1.0]

            P1b = [GCl.RA() - np.cos(GCl.relic_anti_angle)*radius*1.0/np.cos(np.radians(GCl.Dec())), GCl.Dec() + np.sin(GCl.relic_anti_angle)*radius*1.0]
            P2b = [GCl.RA() + np.cos(GCl.relic_anti_angle)*radius*1.0/np.cos(np.radians(GCl.Dec())), GCl.Dec() - np.sin(GCl.relic_anti_angle)*radius*1.0]
            f.show_lines([np.array(zip(P1, P2))], color='w', lw=2., linestyle=':')
            f.show_lines([np.array(zip(P1b, P2b))], color='r', lw=2., linestyle=':')

            if GCl.ratio_relics() > GCl.ratio_relics.vrange[0]:  # Plot if multiple relic
                f.add_label(P1[0], P1[1],  'ratio= %.1e' % (GCl.ratio_relics()), size='x-large', **ciargs)

        # Workaround ... in future it would be better to take the image information from the image and read the contours directly
        try:
            _, center, spixel = maput.FITS2numpy(GCl.mapdic['Raw'])
        except:
            _, center, spixel = 0, (0,0), 7.5

        if relicregions:
            #contours, contourinfos = iom.readDS9regions('Regions/RR_%s.reg'% (GCl.name), spixel, center[0], center[1], pixelcoords=False)
            styles = ['--', ':', '-', '-', '-']
            f.show_polygons( [np.transpose(np.squeeze(region.cnt_WCS)) for region in GCl.regions], lw=2, linestyle=styles[GCl.relics[0].region.rtype.classi+1], **laargs) # , alpha=1.0, facecolor='orange'

        if DS9regions:
            # Load a regions file into APLpy plot
            f.show_regions('Regions/RR_%s.reg' % (GCl.name))

        f.add_scalebar(1)
        f.scalebar.show(1000./kpc, linestyle='solid', linewidth=3., alpha=0.7, **baargs)
        f.scalebar.set_corner('bottom right')
        f.scalebar.set_label('1 Mpc')
        f.scalebar.set_font(size='large', weight='medium', stretch='normal', family='sans-serif',  style='normal', variant='normal')

        if beam:
            f.add_beam()
            f.beam.show()
            f.beam.set_major(GCl.dinfo.beam[0] * u.arcsecond)
            f.beam.set_minor(GCl.dinfo.beam[1] * u.arcsecond)
            f.beam.set_angle(GCl.dinfo.beam[2])  # degrees
            f.beam.set(facecolor=baargs['color'], alpha=0.8, edgecolor='black')
            f.beam.set_color(baargs['color'])

        if show_rot:
            f.add_label(0.97, 0.12, '$\\theta=%.3f$' % (GCl.mockobs.theta), relative=True, style='oblique', size='large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
            f.add_label(0.97, 0.08, '$\\phi  =%.3f$' % (GCl.mockobs.phi)  , relative=True, style='oblique', size='large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)
            f.add_label(0.97, 0.04, '$\\psi  =%.3f$' % (GCl.mockobs.psi)  , relative=True, style='oblique', size='large', horizontalalignment='right', **laargs) #-0.01*len(Cl_name)

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
            f.colorbar.set_axis_label_text(cbar_text)

        """DEVELOPMENT"""
        if vectors:
           pdata = GCl.mapdic['MassSpeed'] #fits with magnitude of signal (use where not enough signal)
           adata = GCl.mapdic['MassAngle'] #fits with angle of signal     (use nan, where no vector should be)
           f.show_vectors(pdata, adata, step=15, scale=1e-2, alpha=0.2, color='blue', lw=2) # , mutation_scale=4 ,ls='-.-', 0.3, head_width=5

#           x  = GCl.mapdic['x'] + 3.0 dsdsd+ GCl.RA  #fits with magnitude of signal (use where not enough signal)
#           y  = GCl.mapdic['y'] + GCl.Dec #fits with angle of signal     (use nan, where no vector should be)
#           dx = GCl.mapdic['dx']  #fits with magnitude of signal (use where not enough signal)
#           dy = GCl.mapdic['dy']  #fits with angle of signal     (use nan, where no vector should be)
#           f.show_arrows(x, y, dx, dy, step=15,scale=1e-2,alpha=0.2, color='blue',lw=2) # , mutation_scale=4 ,ls='-.-', 0.3, head_width=5
        """DEVELOPMENT END"""

        f.set_nan_color((0.5, 0.5, 0.5))

        nowfolder = '%s/Images/' % (survey.outfolder)
        iom.check_mkdir(nowfolder)
        savefile = '%s/%s-%s%s' % (nowfolder, survey.name, GCl.name, 'HR'*highres)
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

    fig, ax = plt.subplots()
    patch = patches.Circle((260, 200), radius=200, transform=ax.transData)
    f._figure.set_clip_path(patch)

    ax.axis('off')
    f.save(savefile+'_circular.png')


def plot_fluxRatio_LAS(surveys):
    """ Test the deviation of the  literatur evalue and measured fluxes of the brightest objects in a cluster """


    for survey in surveys:
        df = survey.fetch_pandas([lambda x: x.GCl.largestLAS, lambda x: x.GCl.flux,
                                        lambda x: x.GCl.flux_lit, lambda x: x.GCl.area],
                                       logs=[False]*3,  keys="dic")
        print('Keys:', df.keys())
        df_clean = df.dropna()
        print(df_clean.keys())

        # use latex for font rendering
        mpl.rcParams['text.usetex'] = True

        plt.clf()
        plt.style.use('default')
        scale = 0.8
        fig, ax = plt.subplots(figsize=(8*scale, 4.7*scale), dpi=200)

        colors = ['b', 'g']
        area = np.power(df_clean['F'], 0.35)*4
        ax.scatter(df_clean['LASmax'], df_clean['F']/df_clean['F_lit'], s=area,  alpha=0.60, color=colors[0], zorder=2)

    farnsw_x = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
                10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5,
                18.0, 18.5, 19.0, 19.5, 20.0]
    farnsw_dec74_pix = [32, 32, 32, 24, 22, 24, 27, 31, 32, 38, 42, 49, 58, 69, 89, 116, 154, 206, 273, 341, 418,
                        494, 570, 641, 706, 766, 820, 868, 912, 952, 986, 1016, 1042, 1066, 1085, 1101, 1114, 1127,
                        1136, 1143, 1148]
    farnsw_dec18_pix = [32, 32, 32, 31, 30, 30, 30, 30, 35, 39, 45, 50, 70, 96, 132, 178, 232, 293, 363, 435, 508,
                        581, 652, 723, 776, 832, 880, 922, 960, 993, 1022, 1047, 1069, 1088, 1106, 1118, 1130, 1140,
                        1147, 1153, 1157]

    farnsw_dec18 = [(1171. - y) / 1140. for y in farnsw_dec18_pix]
    farnsw_dec74 = [(1171. - y) / 1140. for y in farnsw_dec74_pix]

    ax.plot(farnsw_x, [y for y in farnsw_dec18], alpha=0.7, c='grey', zorder=1)
    ax.plot(farnsw_x, [y for y in farnsw_dec74], alpha=0.7, c='grey', zorder=1)
    ax.fill_between(farnsw_x, [y for y in farnsw_dec18], [y for y in farnsw_dec74], color='grey', alpha='0.3', zorder=1)

    powers = [3, 30, 300]
    legl   = [np.power(power, 0.38)*4.5e-0 for power in powers]
    l_iii  = [ax.scatter([],[], s=leg, edgecolors='none', alpha=0.6, color=colors[0]) for leg in legl]
    labels = ['%i' %powe for powe in powers]
    plt.legend(l_iii, labels, ncol=4, frameon=False, fontsize=9, handlelength=1, loc = 1, borderpad=0.4,
               handletextpad=0.2, framealpha=0.70, title='$F_\\mathrm{1.4,\\,NVSS}\\,\mathrm{[mJy]}$', scatterpoints=1)

    ax.set_xlim(0, 20.0)
    ax.set_ylim(ymin=0)
    ax.set_xticks(np.arange(min(ax.get_xlim()), max(ax.get_xlim())+0.5, 3.0))
    ax.set_xlabel('largest $\\mathrm{LAS}\,[\\mathrm{arcmin}]$')
    ax.set_ylabel('$F_\\mathrm{1.4,\\,NVSS} / F_\\mathrm{1.4,\\,lit}$')
    ax.tick_params(direction="in", which='both')
    #ax.set_aspect(1.0/ax.get_data_ratio())

    """
     beam         0.75
     1Mpc z=0.10  8.98
     1Mpc z=0.06 14.28 
     1Mpc z=0.05 16.16
     nominal Largest imagable angular scale by VLA configuration 16.94
    """
    scales = [0.75, 8.98, 16.16, 16.94]
    textl = ['$\\Theta_\\mathrm{FWHM}$','$\\Theta_\\mathrm{z=0.10}$', '$\\Theta_\\mathrm{z=0.05}$', '$\\mathrm{\\Theta_{VLA,D}}$']
    color = ['black', 'b', 'b', 'black']
    height = [ 0.14, 0.2, 0.2, 0.14]
    mod = ((0,0), (0,0), (0,0), (0,0))

    for ii,m in enumerate(scales):
        ax.plot([m]*2, [ax.get_ylim()[0], height[ii] ], '-', c=color[ii], lw=1.8, linestyle=':', alpha=0.7 )
        ax.text(m-0.4+mod[ii][0], height[ii]+0.01+mod[ii][1], textl[ii], fontsize=10, color='black', alpha=0.7)

    #weirdcases = [o.name for o in ClList if (o.flux_lit > 0 and np.log10(o.flux/o.flux_lit) > np.log10(1.3))]
    #print('weirdcases:', weirdcases)

    #fig = plt.figure(figsize=(8 * scale, 4.7 * scale), dpi=200)

    nowfile = 'fluxes_LAS'
    nowfolder = surveys[-1].outfolder + '/PostProcessing/'
    iom.check_mkdir(nowfolder)
    print('Gonna save:  %s' % (nowfolder + nowfile))
    plt.savefig('%s%s.png' % (nowfolder, nowfile), dpi=400)
    plt.savefig('%s%s.pdf' % (nowfolder, nowfile))
    plt.clf()

"""=== Section of single object analysis ==="""
"""========================================="""
def joinpandas(df):
    df_combined = None

    for pdframe in df:
        if df_combined is None:
            df_combined = pdframe
        else:
            df_combined = df_combined.append(pdframe)

    return df_combined


def create_scattermatrix( SurveySamples, plotmeasures, logs=None,  suffix='', shade=True):
    sns.set(style="ticks", color_codes=True)

    """ Creates a scatter matrix, off a list of quantities ... nice! 
    Input: SurveySamples ... there is currently no differnciation between different Survey Samples (symbol-wise or else)
    """

    df = [survey.fetch_pandas(plotmeasures, logs=logs) for survey in SurveySamples]
    original_keys = df[0].keys()
    df_combined = joinpandas(df)
    NSurveys = len(SurveySamples)

    """ Examples of additional plots 
    def hexbin(x, y, color, **kwargs):
         cmap = sns.light_palette(color, as_cmap=True)
         plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)
        
    def corplot(x, y, color, **kwargs):
         cmap = sns.light_palette(color, as_cmap=True)
         plt.imshow(np.abs([x,y].corr()), cmap=cmap, **kwargs) 
    """

    try:
        df_combined['$\\alpha_\mathrm{int}$'] = df_combined['$\\alpha_\mathrm{int}$'].fillna(df_combined['$\\alpha$'])
        df_combined = df_combined.drop(['$\\alpha$'], axis=1)
    except:
        pass

    df_combined = df_combined.reindex(columns=original_keys)
    print('df_combined.Survey.unique()', df_combined.Survey.unique())
    print(df_combined.keys())

    df_combined.to_csv(path_or_buf='/data/Test-%s.csv' % (SurveySamples[0].name))

    g = sns.PairGrid(df_combined, hue="Survey", palette="Set2", dropna=True)
    g = g.map_upper(sns.regplot, scatter_kws={'edgecolors': "white", "linewidth": 1, "alpha": 0.5/np.sqrt(NSurveys)})  #plt.scatter , , edgecolor="white"
    #g = g.map_diag(sns.distplot)
    g = g.map_diag(sns.kdeplot, lw=3, legend=False, alpha=1.0/np.sqrt(NSurveys), shade=True)  #histtype="step"  {'cmap':['Blues_d','Blues']},  ... distplot

    colormaps = ('BuGn', 'Oranges', 'Red') #("Blues", "Blues_d", "Blues_d") #
    #colormaps = sns.cubehelix_palette(8, start=2, rot=0, dark=0, light=.95, reverse=True)

    make_kde.cmap_cycle = cycle(colormaps[0:len(df_combined.Survey.unique())])    #,

    g = g.map_lower(make_kde, alpha=1.0/np.sqrt(NSurveys), shade=shade, shade_lowest=False)

    # from https://stackoverflow.com/questions/52118245/python-seaborn-jointplot-does-not-show-the-correlation-coefficient-and-p-value-o
    for numbered, survey in enumerate(SurveySamples):
        df = survey.fetch_pandas(plotmeasures, logs=logs)
        print('Keys:', df.keys())

        # from https://stackoverflow.com/questions/289.971882/pandas-columns-correlation-with-statistical-significance
        # construct two arrays, one of the correlation and the other of the p-vals
        import pandas as pd
        df_clean = df.dropna()
        rho = df_clean.corr()
        pval = np.zeros([df_clean.shape[1], df_clean.shape[1]])
        for i in range(df_clean.shape[1]):  # rows are the number of rows in the matrix.
            for j in range(df_clean.shape[1]):
                if df_clean.keys()[i] != "Survey" and df_clean.keys()[j] != "Survey":
                    print(i, j, df_clean.shape)
                    JonI = pd.ols(y=df_clean.iloc[:,i], x=df_clean.iloc[:,j], intercept=True)
                    pval[i, j] = JonI.f_stat['p-value']
                    print("Scatterplot, pval for %s-%s" % (df_clean.keys()[i], df_clean.keys()[j]), pval[i, j])

                    """ Only for PhD Thesis, draw line from de gasperin and on fit."""


        xlabels, ylabels = [], []
        for ax in g.axes[-1, :]:
            xlabel = ax.xaxis.get_label_text()
            xlabels.append(xlabel)
        for ax in g.axes[:, 0]:
            ylabel = ax.yaxis.get_label_text()
            ylabels.append(ylabel)

        for i in range(len(xlabels)):
            for j in range(len(ylabels)):
                #g.axes[j, i].xaxis.set_label_text(xlabels[i])
                #g.axes[j, i].yaxis.set_label_text(ylabels[j])
                if i == j:
                    mu, std = norm.fit(df_clean.iloc[:,j])
                    print("Fit results: mu = %.2f,  std = %.2f" % (mu, std))
                #    g.axes[j, i].text(0.5, 0.1+numbered*0.07, '#%i' % (df_clean.shape[0]), zorder=1e10, horizontalalignment='left',
                #                      verticalalignment='center', transform=g.axes[j, i].transAxes)
                if i < j and abs(rho.values[i,j]) > 0.01:
                    g.axes[j, i].text(0.5, 0.1+numbered*0.07, 'correlation: %0.2f' % rho.values[j,i], horizontalalignment='center',
                                      verticalalignment='center', transform=g.axes[j, i].transAxes)
                    slope, intercept, r_value, p_value, std_err = stats.linregress(df_clean.iloc[:,i], df_clean.iloc[:,j])
                    g.axes[j, i].text(0.2, 0.8, "sl=%.2f, ic=%.2f, stde=%.2f" % (slope, intercept, std_err),
                                  horizontalalignment='center',
                                  verticalalignment='center', transform=g.axes[j, i].transAxes)
                    print("Fit results for pairs i,j: sl=%.2f, ic=%.2f, stde=%.2f" % (slope, intercept, std_err))
                    slope, intercept, r_value, p_value, std_err = stats.linregress(df_clean.iloc[:,j], df_clean.iloc[:,i])
                    g.axes[j, i].text(0.2, 0.87, "sl=%.2f,  ic=%.2f, stde=%.2f" % (slope, intercept, std_err),
                                  horizontalalignment='center',
                                  verticalalignment='center', transform=g.axes[j, i].transAxes)
                    print("Fit results for pairs j,i: sl=%.2f,  ic=%.2f, stde=%.2f" % (slope, intercept, std_err))

                if i > j:
                    pass

    #== Save file
    nowfile = 'Scattermatrix'
    nowfolder = SurveySamples[-1].outfolder + '/PostProcessing/'
    iom.check_mkdir(nowfolder)
    print('Gonna save:  %s' % (nowfolder + nowfile))
    plt.savefig('%s%s%s.png' % (nowfolder, nowfile, suffix), dpi=400)
    plt.savefig('%s%s%s.pdf' % (nowfolder, nowfile, suffix))

    plt.clf()


# Taken from https://stackoverflow.com/questions/40726733/plotting-multiple-datasets-on-a-seaborn-pairgrid-as-kdeplots-with-different-colo
def make_kde(*args, **kwargs):
    sns.kdeplot(*args, cmap=next(make_kde.cmap_cycle), **kwargs)


def create_shape_LAS_plot(surveys):
    from scipy.stats import kde
    plt.style.use('default')
    mpl.rcParams['text.usetex'] = True
    plt.rc('text', usetex=True)
    plt.rc('text.latex')

    xmin, xmax = 0.1, 1.42
    ymin, ymax = -1.55, 0

    scale = 0.8
    fig, ax = plt.subplots(figsize=(8 * scale, 4.7 * scale), dpi=200)
    # Create a Rectangle patch

    LAS_line = np.linspace(1, 30, num=50)
    shape_line = np.power(LAS_line, -1.7) * 4.8



    plotmeasures = [lambda x: x.LAS, lambda x: x.iner_rat]


    if 1 == 2:
        df = [survey.fetch_pandas(plotmeasures, logs=[True,True]) for survey in surveys]
        df_combined = joinpandas(df)
        print(df_combined.keys())
        key_LAS = "log$_{10}($LAS [']$)$"
        key_shape = "log$_{10}($$v_\mathrm{PC2}/v_\mathrm{PC1}$$)$"
        data = df_combined[[key_LAS, key_shape]]
        print(data.shape, type(data))
        # Bin sample according to seaborn
        print(data.keys())
        df_NVSS = df_combined[df_combined['Survey'] == 'NVSS']
        df_ELSE = df_combined[df_combined['Survey'] != 'NVSS']


        with sns.axes_style("white"):
            if not df_ELSE.empty:
                sns.jointplot(x=df_ELSE[key_LAS], y=df_ELSE[key_shape], kind="scatter", alpha=0.8, ratio=5);  # color="k",
            g = sns.jointplot(x=df_NVSS[key_LAS], y=df_NVSS[key_shape], kind="scatter", alpha=0.8, ratio=5, color="cornflowerblue")

        # Seaborn figures are square height = X (times X)
        #g.ax_joint.set_xscale('log')
        #g.ax_joint.set_yscale('log')

    if 1 == 1:
        if len(surveys) > 1:
            df = [survey.fetch_pandas(plotmeasures, keys="dic") for survey in surveys[1:]]
            df_combined = joinpandas(df)
            print(df_combined.keys())
            data = df_combined[['LAS', 'iner_rat']]
            print(data.shape, type(data))
            x = data.values[:,0]
            y = data.values[:,1]

            #ax.scatter(x, y, alpha=1 / np.sqrt(len(surveys)), c='salmon')  # , lc='r'

            # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
            nbins = 35
            k = kde.gaussian_kde(data.T)
            xi, yi = np.mgrid[xmin:xmax:nbins * 1j, ymin:ymax:nbins * 1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))

            # contour
            ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.Oranges_r) # plt.cm.BuGn_r, sns.palplot(sns.light_palette("orange", reverse=True))
            ax.contour(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Oranges_r)

            if 1==2:
                g = sns.jointplot(x, y, kind="scatter", color="orange", alpha=0.2, ratio=5)
                g.ax_marg_x.hist(
                    x,
                    alpha=0.5,
                    range=(np.min(x), np.max(x)),
                    color="darkorange",
                    weights=np.ones_like(x) * 1
                )

                g.ax_marg_y.hist(
                    y,
                    orientation='horizontal',
                    alpha=0.5,
                    range=(np.min(y), np.max(y)),
                    color="darkorange",
                    weights=np.ones_like(y) * 1
                )

            # ========

        for survey in [surveys[0]]:
            plotmeasures = [lambda x: x.LAS, lambda x: x.iner_rat]
            df = survey.fetch_pandas(plotmeasures, keys="dic")
            shape = df["iner_rat"]
            LAS = df["LAS"]
            ax.scatter(LAS, shape, alpha=0.8, c='cornflowerblue', zorder=10)  # , lc='r'


    #ax.plot(np.log10(LAS_line), np.log10(shape_line), ls='--', lw=4, alpha=0.5, c="grey")
    ax.tick_params(direction="in", which='both')
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

        #ax.set_xlim([0.2, 1.3])
        #ax.set_ylim([-1.5, 0])
        #ax.set_xticks([2, 3, 5, 7, 10], minor=True)

    #ax.text(0.285, 0.485, "Correlation", transform=plt.gca().transAxes)
    #ax.text(0.56, 0.60, "'Unusual roundish'", transform=plt.gca().transAxes)

    ax.set_xlabel("$\\log_{10}(\mathrm{LAS\,[arcmin])}$")
    ax.set_ylabel("$\\log_{10}(\mathrm{shape}\,s)$")

    nowfile = 'Shape-LAS'
    nowfolder = surveys[-1].outfolder + '/PostProcessing/'
    iom.check_mkdir(nowfolder)
    print('Gonna save:  %s' % (nowfolder + nowfile))
    plt.savefig('%s%s.png' % (nowfolder, nowfile), dpi=400)
    plt.savefig('%s%s.pdf' % (nowfolder, nowfile))
    plt.clf()


def plot_cummulative_flux(surveys, average_relic_count=False):

    plt.style.use('default')
    mpl.rcParams['text.usetex'] = True
    plt.rc('text', usetex=True)
    plt.rc('text.latex')
    scale = 1.0
    fig, ax = plt.subplots(figsize=(6 * scale, 5.5 * scale), dpi=200)
    fig = plt.gcf()
    min_vals, max_vals, cummulatives = [], [], []
    limit = surveys[-1].dinfo.rms*1e3*surveys[-1].relic_filter_kwargs['minrms']

    n_bins = 1200
    bins = np.linspace(np.log10(limit*0.5), np.log10(100000), num=n_bins)
    for survey in [surveys[0]]:
        clusters = survey.FilterCluster(minrel=1)
        fluxes = [np.log10(cl.flux()) for cl in clusters]
        n_relics = [len(cl.filterRelics(**survey.relic_filter_kwargs)) for cl in clusters]

        min_val = min(fluxes)  # min_val = floor(min(data1 + data2))
        max_val = max(fluxes)  # max_val = ceil(max(data1 + data2))
        min_vals.append(min_val)
        max_vals.append(max_val)
        cummulatives.append(len(fluxes))

        # plot the cumulative histogram
        n, bins, patches = ax.hist(fluxes, bins=bins, density=False, histtype='step',
                                   cumulative=-1, color="blue", lw=2, zorder=1000, alpha=0.9)

        if average_relic_count:
            ax2 = ax.twinx()
            ax2.hist(n_relics, bins=bins, density=False, histtype='step',
                     cumulative=-1, color="blue", lw=5, zorder=1000, alpha=0.9, weights=len(clusters))
            ax2.set_ylabel('average relic count', color='darkblue')
            ax2.tick_params('y', colors='darkblue')

    if len(surveys) > 1:
        for survey in surveys[1:]:
            clusters = survey.FilterCluster(minrel=1)
            fluxes = [np.log10(cl.flux()) for cl in clusters]
            n_relics = [len(cl.filterRelics(**survey.relic_filter_kwargs)) for cl in clusters]

            min_val = min(fluxes)  # min_val = floor(min(data1 + data2))
            max_val = max(fluxes)  # max_val = ceil(max(data1 + data2))
            min_vals.append(min_val)
            max_vals.append(max_val)
            cummulatives.append(len(fluxes))

            # plot the cumulative histogram
            n, bins, patches = ax.hist(fluxes, bins=bins,  density=False, histtype='step',
                                       cumulative=-1, color="red", lw=5, alpha=0.2)

        plt.legend([survey.name_short for survey in surveys[0:2]], loc='upper right')

    ax.set_xlim([np.log10(limit), np.max(max_vals) + 0.05])
    ax.set_ylim([0, np.max(cummulatives) + 1])
    ax.set_ylabel('cluster count') #$\sum\mathrm{cluster}\,F_\mathrm{1.4, cluster}>F$
    ax.set_xlabel('$\log_{10}(F\,\mathrm{[mJy]})$')
    plt.tick_params(direction="in", which='both')

    nowfile = 'Flux-distr'
    nowfolder = surveys[-1].outfolder + '/PostProcessing/'
    iom.check_mkdir(nowfolder)
    print('Gonna save:  %s' % (nowfolder + nowfile))
    plt.savefig('%s%s.pdf' % (nowfolder, nowfile))
    plt.clf()


    print('survey.relic_filter_kwargs', survey.relic_filter_kwargs)