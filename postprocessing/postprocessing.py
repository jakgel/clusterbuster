#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 18:21:37 2017

@author: jakobg
"""

    
from __future__ import division, print_function

import os
import corner
import os.path
import numpy as np
import pandas as pd

import clusterbuster.iout.misc as iom
import clusterbuster.iout.surveyplots as ioclass
import clusterbuster.iout.texttables as ioclass_tt
import clusterbuster.dbclasses as dbc
import clusterbuster.surveyclasses as sur
import inference.inference_abcpmc as infer
    
""" Some DEBUGGING and TESTING """
#import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def main():

    create_NVSS = [] #'create_tables', 'create_scattermatrix'] #'create_misc'] #, 'plot_Clusters','correlation template'
    create_ABC_plots = True
    create_PhDplots = False
    create_CWpaperplots = False
    print('###==== PostProcessing: Importing self written .py subroutines ====###')

    colors = ['#0000AA', '#440066', '#AA0000']  # For surveys!!!
    folder = '/data/ClusterBuster-Output/'
    
    print('###==== PostProcessing process ID:', os.getpid())  


    """ load shortlist """
    basepath = os.path.dirname(__file__)
    shortlist = os.path.abspath(os.path.join(basepath, "..", "surveyreal/ClusterList/TableFormReferences.csv"))

    """ In all cases: We load NVSS (for comparison reasons) """
    SurveysSample = []
    NVSSsurvey = iom.unpickleObject('/data/ClusterBuster-Output/%s/pickled/Survey' % ('NVSS'))
    norm = dbc.norm('R200', Nexp=1.0)
    Histo = dbc.Histogram2D(nbins=(64,46), fromto=[[0, 2.*np.pi], [0, 1.5]], norm=norm)
    NVSSsurvey.set_binning(Histo)
    NVSSsurvey.expA = 1.0    #
    NVSSsurvey.expScale = 0.75
    NVSSsurvey.seed_dropout = None
    NVSSsurvey.relic_filter_kwargs = {"Filter": True, "shape":False, "minrms": 8}
    NVSSsurvey.cluster_filter_kwargs = {'minrel': 1, 'zmin': 0.05}


    if 1==2:
        relicsA = NVSSsurvey.fetch_totalRelics()

        # Get alpha and remove nans
        A = np.array([min(-1, relic.alpha()) for relic in relicsA])

        mA = A[~np.isnan(A)].mean()
        print(mA)
        exit()

    if 1==2:
        import copy
        NVSSsurvey_doubleRelics = copy.deepcopy(NVSSsurvey)
        NVSSsurvey_doubleRelics.GCls = [GCl for GCl in NVSSsurvey.GCls if (GCl.gettypestring(vec=True)[2]>0) and (GCl.name != "CIZA J0107")]
        NVSSsurvey = NVSSsurvey_doubleRelics

        for GCl in  NVSSsurvey.GCls:
            if GCl.name == "PLCK G287":
                GCl.M200.value = 18.7e14

        NVSSsurvey.FilterCluster(zmin=0.05) #NVSSsurvey.cluster_filter_kwargs
        for GCl in NVSSsurvey.filteredClusters:
            print("%20s %.2f %6.2f %.2f %.2f" % (GCl.name, np.log10(GCl.P_rest.value), GCl.P_rest.value/1e23, np.log10(GCl.M200.value), GCl.M200.value/1e14))
    if 1==2:
        NVSSsurvey.FilterCluster()
        usethem = NVSSsurvey.filteredClusters
        Nrelics = np.array([len(gcl.filterRelics()) for gcl in usethem])
        M200 = np.array([gcl.M200() for gcl in usethem])
        z = np.array([gcl.z() for gcl in usethem])
        Nrelics = np.expand_dims(Nrelics, axis=1)
        M200 = np.expand_dims(M200, axis=1)
        z = np.expand_dims(z, axis=1)
        measures = np.concatenate((M200, z, Nrelics), axis=1)
        np.save("/data/clustersNVSS_with_mass_redshift", measures)
        exit()

    if 1 == 2:

        print(["%.3f %.3f" % (GCl.R(200)/GCl.R(120), GCl.R(200)/GCl.R(500)) for GCl in NVSSsurvey.GCls])

        NVSSsurvey.GCls = [GCl for GCl in NVSSsurvey.GCls if GCl.name == "A2256"]
        #print(NVSSsurvey.GCls)
        ioclass.plot_Clusters(NVSSsurvey, dynamicscale=False, relicregions=False, DS9regions=False, sectors=False,
                              colorbar=False, infolabel=True, subtracted=True, label_sheme='bright', extralabels=True,
                             filterargs={'zmin': 0.05, 'minimumLAS': 0, 'GClflux': 3.6, 'index': None})
        exit()

    if 1 == 2:
        fetched = NVSSsurvey.fetch_totalRelics()
        print("==== Big. nearby")
        for relic in fetched:
            if relic.Dproj.value/relic.GCl.R200.value < 0.5 and relic.LLS > 1200:
                print(relic.name, relic.Dproj, relic.LLS)
        print("\n ==== Small, distant")
        for relic in fetched:
            if relic.Dproj.value/relic.GCl.R200.value > 0.6 and relic.LLS < 500:
                print(relic.name, relic.Dproj/relic.GCl.R200.value, relic.LLS)
        def cluster_return_polar(Histo, GCl, survey, aligned=False):
            shiftHist = np.roll(Histo.hist.T, -int(aligned * (GCl.relic_pro_index)), axis=1) / survey.AreaHist ** (survey.expA)
            return shiftHist
        exit()

    if 1 == 2:
        if NVSSsurvey.hist_main is not None:
            Histo = NVSSsurvey.hist_main

        NVSSsurvey.set_binning(Histo)
        allHists = None
        for ii, GCl in enumerate(NVSSsurvey.FilterCluster(**NVSSsurvey.cluster_filter_kwargs)):
            """ Deriving the histogram should be a functionality of the survey or the relic cluster, so this should become outdated 
            Beware that at this point, survey.Hist and gcl.Hist are the same objects!
    
            """
            aligned = False
            GCl.updateInformation(Filter=True)
            if GCl.histo is not None and np.sum(GCl.histo.hist) != 0:
                shiftHist = cluster_return_polar(Histo, GCl, NVSSsurvey, aligned=True)
                trueHist  = cluster_return_polar(Histo, GCl, NVSSsurvey, aligned=False)
                shiftHist = np.expand_dims(shiftHist, axis=2)
                trueHist  = np.expand_dims(shiftHist, axis=2)
                if allHists is not None:
                    allHists = np.concatenate((allHists, shiftHist), axis=2)
                    allHists_true = np.concatenate((allHists_true, trueHist), axis=2)
                else:
                    allHists = shiftHist
                    allHists_true = trueHist
        np.save("/data/allHist_NVSS_full_polar_aligned",  allHists)
        np.save("/data/allHist_NVSS_full_polar",  allHists_true)

    if 1 == 2:
        def compute_moments(sparseA, sparseD, sparseW):
            x = np.cos(sparseA) * sparseD
            y = np.sin(sparseA) * sparseD

            moment_00 = np.sum(sparseW)
            moment_10 = np.sum(sparseW * x)
            moment_01 = np.sum(sparseW * y)
            moment_11 = np.sum(sparseW * x * y)
            moment_20 = np.sum(sparseW * x * x)
            moment_02 = np.sum(sparseW * y * y)


            _xm = moment_10 / moment_00
            _ym = moment_01 / moment_00

            moment_central_11 = moment_11 / moment_00 - _xm * _ym
            moment_central_20 = moment_20 / moment_00 - _xm * _xm
            moment_central_02 = moment_02 / moment_00 - _ym * _ym
            moment_angle = .5*np.arctan2(2*moment_central_11, moment_central_20-moment_central_02)
            return (moment_angle, moment_00, moment_10, moment_01, moment_11  - _xm * _ym , moment_20 - _xm * _xm, moment_02 -_ym * _ym)

        arrays=[]
        angles = [ii/100*2*np.pi for ii in range(20)]
        for angle in angles:
            sparseA = [angle]*2
            sparseD = [1,0.9]
            sparseW = [1]*2
            array = compute_moments(sparseA, sparseD, sparseW)
            arrays.append(array)

    if 1 == 2:
        filtered_relics = NVSSsurvey.fetch_totalRelics()
        flux = np.asarray([relic.flux() for relic in filtered_relics])
        alpha_error = np.asarray([relic.alpha.std[0] for relic in filtered_relics])
        alpha_error[alpha_error == 0.3] = np.nan
        np.save("/data/NVSS_flux", flux)
        np.save("/data/NVSS_alpha_error", alpha_error)

    if 1 == 2:
        GCls = []
        import clusterbuster.surveyclasses as sc
        import clusterbuster.sourceextraction as se

        angles_put = [ii/20*2*np.pi for ii in range(20)]

        for ii, angle_test in enumerate(angles_put):
            print('___angle:', angle_test*180/np.pi)
            dist = 15*(1+ii/9)
            y_pos = np.cos(angle_test)*dist
            x_pos = np.sin(angle_test)*dist


            detinfo = sc.DetInfo(name='test', beam=[7.5, 7.5, 0], spixel=7.5, rms=0, limit=0.1, nucen=1.4, center=[0,0], pcenter=[50, 50], survey='NVSS', telescope='VLA-D')
            rinfo = sc.RelicRegion(name='', cnt=[], rtype=1, alpha=1.2)
            xx, yy = np.mgrid[:100, :100]
            ellipse = 1*(xx - 50 - x_pos) ** 2 + 1*(yy - 50 - y_pos) ** 2
            donut = (ellipse < 30) #& (ellipse >= 0)
            image = donut*10.
            relics = se.RelicExtraction(image, 0.2, dinfo=detinfo, rinfo=rinfo)
            GCl = sc.Galaxycluster("%i" % ii, 0, 0, 0.2, relics=relics)
            GCls.append(GCl)

    # just because there is not a single full-hist-option ... or cluster hist option
    if 1 == 2:
        Histo = NVSSsurvey.hist_main
        NVSSsurvey.set_binning(Histo)
        NVSSsurvey.FilterCluster(**NVSSsurvey.cluster_filter_kwargs)

        expA = 1
        allHists = None
        for gcl in NVSSsurvey.filteredClusters:
            shiftHist = cluster_return_polar(Histo, gcl, NVSSsurvey)
            shiftHist = np.expand_dims(shiftHist, axis=2)
            if allHists is not None:
                allHists = np.concatenate((allHists, shiftHist), axis=2)
            else:
                allHists = shiftHist
        print(allHists.shape, np.sum(allHists.flatten()))
        np.save("/data/allHist_NVSS",  allHists)


        _, comprad, comppol, _, _ = NVSSsurvey.polar(zmin=0.05, conv=0)
        np.save("/data/comprad", comprad)
        np.save("/data/compol", comppol)


    def plot_correlations(data, suffix):
        from string import ascii_letters
        import pandas as pd
        import seaborn as sns
        import copy

        sns.set(style="white")

        data_frame = pd.DataFrame(data=data,  # values
                                  columns=labelset)  # 1st row as the column names

        # Compute the correlation matrix
        corr = data_frame.corr()

        # Generate a mask for the upper triangle
        mask = np.zeros_like(corr, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True

        # Set up the matplotlib figure
        f, ax = plt.subplots(figsize=(11 * 0.7, 9 * 0.7))

        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(220, 10, as_cmap=True)

        # Draw the heatmap with the mask and correct aspect
        #corr_rounded = copy.deepcopy(corr)
        #corr_rounded.iloc[:, :] = np.round(corr_rounded.values, decimals=2)


        sns.heatmap(corr, mask=mask, cmap=cmap, vmin=.65, vmax=.65, center=0,
                    square=True, linewidths=.5, annot=True, fmt="+.2f", annot_kws={'size':12*8/(corr.values.shape[0]+2)},
                    cbar_kws={"shrink": .6, "orientation": "horizontal", "label": "correlation"})

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")
        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_yticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        # for i in range(len(data_frame.keys())):
        #    for j in range(len(data_frame.keys())):
        #        if i>j:
        #            text = ax.text(j, i, "%+.2f" % corr.values[i, j],
        #                           ha="center", va="center", color="w")

        plt.savefig('%s/correlation-matrix-%s.png' % (folderN, suffix))
        plt.savefig('%s/correlation-matrix-%s.pdf' % (folderN, suffix))

    if create_NVSS:
        """ NVSS """
        #        ioclass.plot_Clusters(NVSSsurvey, dynamicscale=False, relicregions=True, DS9regions=False, sectors=False, colorbar=False, subtracted=False, shapes=False)
        if 'create_misc':
            #ioclass.plot_cummulative_flux([NVSSsurvey])
            ioclass.plot_fluxRatio_LAS([NVSSsurvey])
            #ioclass.create_shape_LAS_plot([NVSSsurvey])

        if "correlation template" in create_NVSS:
            mpl.rcParams['text.usetex'] = True
            plt.style.use('default')

            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            scale = 0.8
            fig, ax = plt.subplots(figsize=(8 * scale, 4.7 * scale), dpi=200)

            for gcl in NVSSsurvey.FilterCluster(**NVSSsurvey.cluster_filter_kwargs):
                gcl.updateInformation(Filter=True)
                gcl.relics_polarDistribution()
                arrays = [val / max(gcl.fitarray_pro) for val in gcl.fitarray_pro]
                ticks = [tick for tick in gcl.histo.ticks[0]]
                ticks_relative = [((((tick - gcl.histo.ticks[0][
                    gcl.relic_pro_index]) + 1 / 2 * np.pi) % np.pi) - 1 / 2 * np.pi) * 180 / np.pi for tick in
                                  gcl.histo.ticks[0]]

                ticks_rolled, arrays_rolled = zip(*sorted(zip(ticks_relative, arrays)))
                ax.plot(ticks_rolled, arrays_rolled, alpha=0.4, c='cornflowerblue', lw=1.5)  # , lc='r'
                # ax.scatter(ticks_rolled, arrays_rolled, alpha=0.1, c='r')  # , lc='r'

                ticks_example = np.linspace(-90, 90, 100)
                ax.plot(ticks_example, np.abs(np.cos(ticks_example / 180 * np.pi)), c='black')

            ax.tick_params(direction="in", which='both')
            ax.set_xlim([-90, 90])
            ax.set_ylim([0, 1])
            plt.xlabel("$\phi_\\mathrm{corr}-\phi_\\mathrm{best}\;[^\circ]$")
            plt.ylabel("$R(\phi_\\mathrm{corr}) / R(\phi_\\mathrm{best})$")
            plt.savefig("/data/ClusterBuster-Output/NVSS/PostProcessing/Correlation-template.pdf")

        if 'create_tables' in create_NVSS:
            ioclass_tt.GClList2table_paper('/home/jakobg/Dropbox/PhDThesis/tables/GClTable', NVSSsurvey, longtab=False, shortlist=shortlist)
            ioclass_tt.RList2table_paper('/home/jakobg/Dropbox/PhDThesis/tables/RelicsTable', NVSSsurvey, longtab=False)
            #ioclass.plot_RelicEmission_polar(NVSSsurvey, additive=True, single=True, mirrored=False, plottype='flux',
            #                                 cbar=True, aligned=False, dpi=600, title=None, add_pi=0.5)

        if 'create_scattermatrix' in create_NVSS:
            newdist = lambda x: dbc.measurand(x.Dproj_pix() / x.GCl.R200(), 'Dproj', label='$D_\mathrm{proj,rel}$',
                                              un='$R_{200}$')
            newalpha = lambda x: dbc.measurand(np.abs(x.alpha()), '|alpha|', label='$|\\alpha|$', un=None)
            logs_alpha = [True, True, False, True, False, True]
            logs_alpha_new  = [True, True, False, False]
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.P_rest, lambda x: x.M200, lambda x: x.z], logs=[True,True,False], suffix='_z_Gcls', gcls=True)
            exit()
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.alpha, lambda x: x.Dproj_pix], logs=[False,True])
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.alpha, newdist], logs=[False,True], suffix='_PhDplot' )
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.P_rest, lambda x: x.GCl.M200], logs=[True,True], suffix='_PhDplot')
            #ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.GCl.z, lambda x: x.Mach], suffix='_Mach-z')
            #ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, lambda x: x.Dproj_pix], logs=logs_alpha, suffix='_large')
            #ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.Mach, lambda x: x.Dproj_pix], suffix='_large_mach')
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist], logs=logs_alpha, suffix='_newdist')
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist], logs=logs_alpha_new, suffix='_newdistB')
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist, lambda x: x.GCl.z, lambda x: x.GCl.M200], logs=logs_alpha, suffix='_OWH')
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, lambda x: x.Dproj_pix, lambda x: x.GCl.z, lambda x: x.GCl.M200], logs=logs_alpha, suffix='_OWH2')
            ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.GCl.z], logs=logs_alpha, suffix='_z')
            #ioclass.create_scattermatrix([NVSSsurvey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.GCl.z, lambda x: x.GCl.M200], logs=logs_alpha, suffix='_M200')

        #print('len(NVSSsurvey.FilterCluster(minrel=1))', len(NVSSsurvey.FilterCluster(minrel=1))
        if 'plot_Clusters' in create_NVSS:
            ioclass.plot_Clusters(NVSSsurvey, dynamicscale=False, relicregions=False, DS9regions=False, sectors=False, colorbar=False, infolabel=True,
                                  subtracted=True, filterargs={'zmin': 0.05, 'minimumLAS': 0, 'GClflux': 3.6, 'index': None})




    def add_correlation(corner_plot, x, text=('top')):
        import pandas as pd
        df_clean = pd.DataFrame(x)
        rho = df_clean.corr()
        rows = rho.shape[0]
        for i in range(rows*rows):
            if i % rows < int(i/rows):
                corner_plot.axes[i].text(0.05, 0.05, '%+.2f' % (rho.values[i % rows, int(i/rows)]),
                                         transform=corner_plot.axes[i].transAxes,
                                         size=16) #, color='r'


    def add_literature_values(corner_plot, id_B0=1, id_kappa=2):
        df_clean = pd.DataFrame(x)
        values = [ (0,0,5), (0.5,1.0), (0.3,0.7)]
        rows = rho.shape[0]
        for i in range(rows*rows):
            if i % rows < int(i/rows):
                corner_plot.axes[i].text(0.05, 0.05, '%+.2f' % (rho.values[i % rows, int(i/rows)]),
                                         transform=corner_plot.axes[i].transAxes,
                                         size=16) #, color='r'

    if create_ABC_plots:


        #plt.style.use('classic')
        allpools  = []
        plottings = []
        model_samples_total = 0
        mode = "thesis_largestRun" #"4vs4_new" #
        folderN = '/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_Run_16'
        show_metrics = False
        metric_log = False


        if mode == "1vs1":
            labelset = ['$\\log_{10} \Delta_\\mathrm{count}$', "$\log_{10}(\\xi_\\mathrm{e}$)"]
            labelset = ['$\Delta_\\mathrm{count}$', "$\log_{10}(\\xi_\\mathrm{e}$)"]
            rangeset = [[-0.03,  0.3],
                        [-4.9, -4.4]
                        ]
        elif mode == "effic-B-tiny":

            labelset = [r'$\log_{10} \Delta_\mathrm{count}$', r'$\log_{10} \Delta_\mathrm{average\,relic}$',
                        r"$\log_{10} \xi_e$", r"$log_{10} B_0$"]
            rangeset = [[-1, 1],
                        [-6, -3.5],
                        ]
        elif mode == "3vs3_old":
            labelset = [r'$\log_{10} \Delta_\mathrm{count}$', r'$\log_{10} \Delta_\mathrm{average\,relic}$',
                        r'$\log_{10} \Delta_\mathrm{2DKS}$', r"$\log_{10} \xi_e$", r"$log_{10} B_0$", r"$\kappa$"]
            rangeset = [[-0.3, 0.0],
                        [-1.7, 0.4],
                        [-0.6, 0.5],
                        [-6, -3.5],
                        [-0.8, 2.0],
                        [-1.2, 1.0]
                        ]
        elif mode == "3vs3_new":
            labelset = [r"$\log_{10} \xi_e$", r"$\log_{10} B_0$", r"$\kappa$"]
            rangeset = [[-6, -4.0],
                        [-0.5, 3.0],
                        [-0.0, 3.0]
                        ]
        elif mode == "4vs3_new":
            labelset = ["$\Delta_\mathrm{count}$", "$\Delta_\mathrm{average\,relic}$", "$\Delta_\mathrm{2DKS}$",
                        r"$\log_{10} \xi_e$", r"$\log_{10} B_0$", r"$\kappa$", "$b_\mathrm{pca\,filter}$"]
            rangeset = [[0, 2.0],
                        [0, 1.0],
                        [0, 1.0],
                        [-5.5, -4.0],
                        [-0.5, 2.5],
                        [-0.5, 2.0],
                        [-0.3,0.3]
                        ]
            #eps_use_wide = [4 ,2, 1]
            eps_use_wide = [1.39420435, 0.41220674, 0.90492575]
            #eps_use_narrow = [0.35668561, 0.31255057, 0.74173276]
            eps_use_narrow = [0.62275237, 0.29084472, 0.695507]  # Last iteration
            eps_select = lambda x: [float(x[0]), float(x[1]), float(x[2]),
                                     np.log10(float(x[4])), np.log10(float(x[5])), float(x[6]), float(x[8])]
            #par_select = lambda x: (float(x[5]) < 1.0 and float(x[6]) > 0.0) and float(x[8]) > -0.3
            par_select = lambda x: float(x[5]) < 5 #and float(x[8]) > -0.0
            #par_select = lambda x: float(x[8]) > -0.0 #float(x[8]) < 0.10 and
        elif mode == "4vs4_new":
            labelset = ["$\Delta_\mathrm{count}$", "$\Delta_\mathrm{average\,relic}$", "$\Delta_\mathrm{2DKS}$", "$\Delta_\mathrm{alpha,mean}$",
                        r"$\log_{10} \xi_e$", r"$\log_{10} B_0$", r"$\kappa$", "$b_\mathrm{pca\,filter}$"]
            rangeset = [[0, 2.0],
                        [0, 1.0],
                        [0, 1.0],
                        [0, 0.1],
                        [-5.5, -4.0],
                        [-0.5, 2.5],
                        [-0.5, 2.0],
                        [-0.5,0.5]
                        ]
            eps_use_wide = [3, 2, 1, 0.3]
            eps_use_narrow = [ 0.7989435,  0.26666418,  0.71958977,  0.01344388]  #0.35, 0.38, 0.62, 0.03]
            #eps_use_narrow = [0.35668561, 0.31255057, 0.74173276]
            eps_select = lambda x: [float(x[0]), float(x[1]), float(x[2]), float(x[3]),
                                     np.log10(float(x[5])), np.log10(float(x[6])), float(x[7]), float(x[9])]
            par_select = lambda x: float(x[6]) < 5
        elif mode == "3vs5_new":
            labelset = ["$\Delta_\mathrm{count}$", "$\Delta_\mathrm{average\,relic}$", "$\Delta_\mathrm{alpha,mean}$",
                        r"$\log_{10} \xi_e$", r"$\log_{10} B_0 [\mu G]$", r"$\kappa$",
                        r"$\log_{10}\,\Xi_\mathrm{ratio}$", r"$\log_{10} t_1 \mathrm{[Gyr]}$", r"$\log_{10} t_2 \mathrm{[Gyr]}$"]
            rangeset = [[0, 2.0],
                        [0, 1.0],
                        [0, 0.1],
                        [-6, -4.0],
                        [-0.8, 3.0],
                        [-0.5, 3.0],
                        [-5, 0.5],
                        [-1.3, 2.0],
                        [-1.3, 2.0]
                        ]

            eps_use_wide = [0.95, 4, 0.1]
            eps_use_narrow = [0.4, 0.5, 0.02]
            eps_select = lambda x: [float(x[0]), float(x[1]), float(x[2]),
                                     np.log10(float(x[4])), np.log10(float(x[5])), float(x[6]),
                                     np.log10(float(x[7])), np.log10(float(x[8])), np.log10(float(x[9]))]

        elif mode == "thesis_largestRun":
            labelset = ["$\Delta_\mathrm{count}$", "$\Delta_\mathrm{average\,relic}$", "$\Delta_\mathrm{2DKS}$", "$\Delta_\mathrm{alpha,mean}$",
                        r"$\log_{10} \xi_e$", r"$\log_{10} B_0 [\mu G]$", r"$\kappa$", "$b_\mathrm{pca\,filter}$",
                        r"$\log_{10}\,\Xi_\mathrm{ratio}$", r"$\log_{10} t_1 \mathrm{[Gyr]}$", r"$\log_{10} t_2 \mathrm{[Gyr]}$"]

            rangeset = [[0, 2.0],
                        [0, 1.0],
                        [0, 1.0],
                        [0, 0.1],
                        [-5.5, -4.0],
                        [-0.5, 2.5],
                        [-0.5, 2.0],
                        [-0.3, 0.3],
                        [-5, 0.5],
                        [-1.3, 1.0],
                        [-1.3, 1.0]
                        ]

            eps_use_wide = [3, 0.7, 0.98, 0.1]
            eps_use_narrow = [1.13, 0.30, 0.835, 0.01966]
            eps_use_narrow = [0.98398749, 0.27925968, 0.78109853, 0.01610189]
            #eps_use_narrow = [0.63, 0.25, 0.78, 0.01]
            eps_select = lambda x: [float(x[0]), float(x[1]), float(x[2]), float(x[3]),
                                     np.log10(float(x[5])), np.log10(float(x[6])), float(x[7]), float(x[14]),
                                     np.log10(float(x[10])), np.log10(float(x[8])), np.log10(float(x[9]))]
            par_select = lambda x: float(x[14]) > -0.0
        if not show_metrics:
            labelset = labelset[len(eps_use_wide):]
            rangeset = rangeset[len(eps_use_wide):]



        quantiles = [0.16, 0.5, 0.84]
        levels = [0.90, 0.60, 0.30] #None #quantiles #(1 - np.exp(-0.16), 1 - np.exp(-0.5), 1 - np.exp(-0.84),)
        used_cmap = plt.get_cmap('viridis')    # hot, viridis, magma, hot, 'viridis_r'
        kwargs = {"quantiles": quantiles,  "labels": labelset, "show_titles": True,
                  "label_kwargs": {"fontsize": 18}, "title_kwargs": {"fontsize": 18},
                  "contourf_kwargs": {'cmap': used_cmap, 'colors': None},  "range": rangeset,
                  "hist2d_kwargs": {"levels": levels}}
        kwargs2 = {"quantiles": quantiles, "levels": levels, "labels": labelset, "show_titles": True,
                   "label_kwargs": {"fontsize": 18}, "contourf_kwargs": {'cmap': used_cmap, 'colors': None},
                   "title_kwargs":{"fontsize": 18, "loc":'left'}, "hist2d_kwargs": {"levels": levels}}

        if 1==1:
            all_selected = []
            for line in open(folderN + '/logfile.txt'):
                line = line.rstrip('\n')
                listWords = line.split(" ")
                listWords = [i for i in listWords if i]
                all_selected.append(eps_select(listWords))
            all_selected = np.asarray(all_selected)
            np.save('%s/Converted_logfile.npy' % (folderN), all_selected)


        launch_pools = iom.unpickleObjectS(folderN+'/launch_pools.pickle')
        for pools in launch_pools:
            infer.plot_abctraces(pools, folderN)
            print('len(pools):', len(pools))

            for ii, pool in enumerate(pools):
                print('create_ABC_plots:pool.dists.shape', pool.dists.shape)
                print('create_ABC_plots:pool.thetas.shape', pool.thetas.shape)
                print('create_ABC_plots:pool.eps', pool.eps)
                print('create_ABC_plots:pool.ratio', pool.ratio)
                print('create_ABC_plots:pool.dists[0:2]', pool.dists[0:2])
                print('create_ABC_plots:pool.thetas[0:2]', pool.thetas[0:2])
                print('create_ABC_plots:pool.ws[0:2]', pool.ws[0:2])
                print('create_ABC_plots:pool.dists[-1]', pool.dists[-1])
                print('create_ABC_plots:pool.thetas[-1]', pool.thetas[-1])
                model_samples_total += int(pool.dists.shape[0]/pool.ratio)

                if show_metrics:
                    if metric_log:
                        plotting_old = np.concatenate((np.log10(pool.dists), pool.thetas), axis=1)
                    else:
                        plotting_old = np.concatenate((pool.dists, pool.thetas), axis=1)
                else:
                    plotting_old = pool.thetas

                print(plotting_old.shape)
                allpools.append(pool)
                plottings.append(plotting_old)
        print('len(pools)', len(pools))


        print('model_samples_total:', model_samples_total)
        print(np.asarray(plottings).shape)


        #================
        #import triangle
        import abcpmc

        if 1==2:
            samples = None
            for pool in pools:

                print(pool.eps, pool.dists.shape, pools[-1].eps)
                selected = np.ones_like(pool.thetas[:,0])
                print(pool.dists.shape, selected.shape)
                for ii in range(pool.dists.shape[1]):
                    one_criterion= pool.dists[:,ii] < pools[-1].eps[ii]
                    selected = selected * one_criterion
                print(selected.shape)
                selected = np.where(selected==1)
                print(type(selected), pool.thetas.shape)
                selected_thetas = pool.thetas[selected]
                selected_weight = pool.ws[selected]

                if samples is None:
                    samples = selected_thetas
                    weights = selected_weight
                else:
                    samples = np.concatenate((samples, selected_thetas), axis=0)
                    weights = np.concatenate((weights, selected_weight), axis=0)

            print(samples.shape, weights.shape)
            fig = corner.corner(samples, weights=weights, **kwargs2)
            add_correlation(fig, samples)
            plt.savefig("/data/Test.pdf")


            for mean, std in zip(*abcpmc.weighted_avg_and_std(samples, weights, axis=0)):
                print(u"mean: {0:>.4f} \u00B1 {1:>.4f}".format(mean, std))

            exit()





        print('eps_use_wide', eps_use_wide)
        all_selected = []
        for line in open(folderN + '/logfile.txt'):
            line = line.rstrip('\n')
            listWords = line.split(" ")
            listWords = [i for i in listWords if i]

            if all(float(listWords[ii]) < eps_use for ii, eps_use in enumerate(eps_use_wide)):
                all_selected.append(eps_select(listWords))
        all_selected = np.asarray(all_selected)
        print(all_selected.shape)
        print(all_selected[:, 0].flatten().shape)
        for ii in range(len(eps_use_wide)):
            all_selected_hist = all_selected[:,ii].flatten()
            plt.close("all")
            plt.hist(all_selected_hist, cumulative=False, bins=40)
            plt.savefig('%s/Metric_Histo_%i.pdf' % (folderN, ii))
            plt.close("all")
            plt.hist(all_selected_hist, cumulative=True, bins=3000)
            plt.savefig('%s/Metric_CumHisto_%i.pdf' % (folderN, ii))
            plt.close("all")

        # WideEpsilon
        if not show_metrics:
            all_selected = all_selected[:, len(eps_use_wide):]
        plottet_corner = corner.corner(all_selected, **kwargs2)
        add_correlation(plottet_corner, all_selected)
        plt.savefig('%s/WideEpsilon.pdf' % (folderN))
        plt.savefig('%s/WideEpsilon.png' % (folderN))
        plt.close("all")

        all_selected = []
        for line in open(folderN + '/logfile.txt'):
            line = line.rstrip('\n')
            listWords = line.split(" ")
            listWords = [i for i in listWords if i]
            if all(float(listWords[ii]) < eps_use for ii,eps_use in enumerate(eps_use_narrow)):
                all_selected.append(eps_select(listWords))
        all_selected = np.asarray(all_selected)
        if not show_metrics:
            all_selected = all_selected[:, len(eps_use_wide):]
        print(all_selected[:, 0].flatten().shape)

        # SecondLastEps
        plottet_corner = corner.corner(all_selected, **kwargs2)
        add_correlation(plottet_corner, all_selected)
        plt.savefig('%s/SecondLastEps.pdf' % (folderN))
        plt.savefig('%s/SecondLastEps.png' % (folderN))
        plt.close("all")



        all_selected_priors = []
        for line in open(folderN + '/logfile.txt'):
            line = line.rstrip('\n')
            listWords = line.split(" ")
            listWords = [i for i in listWords if i]
            if all(float(listWords[ii]) < eps_use for ii, eps_use in enumerate(eps_use_narrow)):
                if par_select(listWords):
                    all_selected_priors.append(eps_select(listWords))
        all_selected_priors = np.asarray(all_selected_priors)
        if not show_metrics:
            all_selected_priors = all_selected_priors[:,len(eps_use_wide):]
        print(all_selected_priors[:, 0].flatten().shape)

        # SecondLastEps
        plottet_corner = corner.corner(all_selected_priors, **kwargs2)
        add_correlation(plottet_corner, all_selected_priors)
        plt.savefig('%s/SecondLastEps_with_priors.pdf' % (folderN))
        plt.savefig('%s/SecondLastEps_with_priors.png' % (folderN))
        plt.close("all")


        for ii, plotting in enumerate(plottings):
            print(np.asarray(plotting).shape)
            corner.corner(np.asarray(plotting), **kwargs)
            plt.savefig('%s%i.pdf' % (folderN, ii))
            plt.savefig('%s%i.png' % (folderN, ii))

        # Ratio development
        ratios = np.asarray([pool.ratio for pool in pools])
        print(ratios.shape)
        plt.plot(ratios[1:]/ratios[0:-1])
        plt.xlabel('Iteration')
        plt.ylabel('Relative change of acceptance rate since last iteration')
        plt.savefig('%s/Acceptance_development.pdf' % (folderN))
        plt.savefig('%s/Acceptance_development.png' % (folderN))
        plt.close("all")

        #LastThree
        if len(pools) >= 3:
            corner.corner(np.concatenate((plottings[-3], plottings[-2], plottings[-1]), axis=0), **kwargs)
            plt.savefig('%s/LastThree.pdf' % (folderN))
            plt.savefig('%s/LastThree.png' % (folderN))

        #LastTwo
        if len(pools) >= 2:
            corner.corner(np.concatenate((plottings[-2], plottings[-1]), axis=0), **kwargs)
            plt.savefig('%s/LastTwo.pdf' % (folderN))
            plt.savefig('%s/LastTwo.png' % (folderN))
            plt.close("all")

        # Generate the pandas dataframe
        for hh, data in enumerate(plottings):
            plot_correlations(data, '%02i' % (hh))
        plot_correlations(all_selected, "SecondLastEps")
        plot_correlations(all_selected_priors, "SecondLastEps_priors")


    if create_PhDplots:
        """ MUSIC-2 """
        surveynames = []
        mainname    = "Survey"  #'Threehundrets'
        #========== My specific cohord plot =============#
        plt.style.use('ggplot')

        proclist = range(0, 36)  # range(0,36), 313
        #proclist = range(5600, 5600+50)  # range(1360, 3288+250)  '/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_Run_11/'
        folder = '/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_TestingSimpleParas_Nuza_woPhaseFilter_0.75_100kpc_10times_sensitivity/'
                #'/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_TestingSimpleParas_Nuza_woPhaseFilter_0.75_20kpc/'
                #'/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_Run_13/'
                #
                # '/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_TestingSimpleParas_Nuza_woPhaseFilter_0.75_100kpc_3times_sensitivity/surveys/'
                # '/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_TestingSimpleParas_Nuza_woPhaseFilter_0.75_100kpc_10times_sensitivity/surveys/'
                 #'/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_TestingSimpleParas_Nuza_woPhaseFilter_0.75_100kpc/surveys/'
        from os.path import normpath, basename
        outname = basename(normpath(folder))
        folder += "surveys/"

        def radialstatistics(surveypaths,positive=True, normalize=True, density=True,color='red'):
            fullarray = []
            
            for surveypath in surveypaths:
                try:
                    survey = iom.unpickleObject(surveypath)
                    print(survey.name)
                except:
                    continue
                survey.relic_filter_kwargs = {"Filter": True, 'minrms': 8, "shape_pca": True}
                survey.dinfo = survey.GCls[0].dinfo       
                survey.scatterkwargs = {"alpha": 0.15, "fmt": "^", "markersize": 7}
                survey.histkwargs = {"alpha": 0.15}
                survey.set_seed_dropout()
                survey.set_binning(Histo)


                (radial, radial_tickz) = survey.polar(positive=positive, normalize=normalize, density=density)[1]
                plt.plot(radial_tickz, radial,  color=color, alpha=0.02)
                plt.fill_between(radial_tickz, radial, where=radial>=0, color=color, alpha=0.008)
                fullarray.append(radial)
                
            fullarray = np.stack([np.asarray(i) for i in fullarray], axis=0)
            plt.plot(radial_tickz, np.median(fullarray, axis=0), color=color, alpha=0.8)

           
        def CombinedSurvey_Statistics():
            surveypaths = []
            plt.style.use('default')
            scaling = 1
            fig, (ax) = plt.subplots(1, 1, figsize=(7 * scaling, 4.2 * scaling), dpi=300)
            for ii in proclist:
                    surveypaths.append('%s/Survey_%05i' % (folder, ii))


            for positive in [False,True]:
                density = True
                normalize = True
                
                radialstatistics(surveypaths, positive=positive)
                (radial,radial_tickz) = NVSSsurvey.polar(positive=positive, normalize=normalize, density=density)[1]

                ax.plot(radial_tickz, radial, color='blue', alpha=0.3)
                ax.fill_between(radial_tickz, radial, where=radial>=0, color='blue', alpha=0.2)

                patch1 = mpatches.Patch(color='blue', label='NVSS', alpha=0.2)
                patch2 = mpatches.Patch(color='red', label='synthetic survey', alpha=0.2)
                plt.legend(handles=[patch1, patch2])

                if positive:
                    plt.xlabel('Distance to cluster centre $[R_{200}]$')
                    ax.set_xlim(left=0, right=1.5) #ymax=0.05
                else:
                    plt.xlabel('Distance along pro-relic axis $[R_{200}]$')
                    ax.set_xlim(left=-1.5, right=1.5) #ymax=0.05
                plt.ylabel('Signal density $S\,[R_{200}^{-1}]$')
                ax.set_ylim(bottom=0) #ymax=0.05

                ax.tick_params(direction="in", which='both')
                ax.tick_params(direction="in", which='major', right=False, top=True, labelright=False)

                plt.savefig('/data/ClusterBuster-Output/Signal_distance%s.pdf' % ('_radial' * positive))
                plt.clf()
            return 0

        for i in proclist:
                surveyname = '%s_%05i' % (mainname,i)
                surveynames.append(surveyname)

        if 1 == 2:
            CombinedSurvey_Statistics()
            exit()

        if 1 == 2:
            total = None
            for surveyname in surveynames:
                print('%s/%s' % (folder, surveyname))
                try:
                    survey = iom.unpickleObject('%s/%s' % (folder, surveyname))
                except:
                    continue
                print(survey.name)
                survey.relic_filter_kwargs = {"Filter": True, 'minrms': 8, "shape_pca": True}
                survey.set_seed_dropout()
                relicsA = survey.fetch_totalRelics()

                # Get alpha and remove nans
                A = np.array([ [relic.P_rest.value, relic.LLS.value] for relic in relicsA])
                if total is None:
                    total = A
                else:
                    total = np.concatenate((total,A), axis=0)
                    print(total.shape)
            np.save("/data/SurveyLLSPower_10times", total)
            exit()

        SurveysSample_full = [NVSSsurvey]
        if 1 == 2:
            for surveyname in surveynames:
                print('%s/%s' % (folder, surveyname))
                try:
                    survey = iom.unpickleObject('%s/%s' % (folder, surveyname))
                except:
                    continue
                print(survey.name)

                survey.set_seed_dropout()
                survey.dinfo = survey.GCls[0].dinfo
                survey.scatterkwargs = {"alpha": 0.15, "fmt": "^", "markersize": 7}
                survey.histkwargs = {"alpha": 0.15}
                survey.set_binning(Histo)
                survey.relic_filter_kwargs.update({"Filter": True, 'minrms': 8, "shape_pca": True})

                SurveysSample_full.append(survey)


                newdist = lambda x: dbc.measurand(x.Dproj_pix() / x.GCl.R200(), 'Dproj', label='$D_\mathrm{proj,rel}$',
                                                  un='$R_{200}$')
                plotmeasures = [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist, lambda x: x.GCl.z,
                                lambda x: x.GCl.M200, lambda x: x.flux, lambda x: x.shape_advanced(), lambda x: x.LAS,
                                lambda x: x.iner_rat]
                logs = [True, True, False, True, False, True, True, False, True, True]

                #for survey in SurveysSample_full:
                #    survey.relic_filter_kwargs = {"Filter":True, "shape":False, "minrms":8}
                #df = [survey.fetch_pandas(plotmeasures, logs=logs) for survey in SurveysSample_full]
                #original_keys = df[0].keys()
                #df_combined = ioclass.joinpandas(df)
                #df_combined.to_csv("/data/surveys_compressed.csv")

            ioclass.create_shape_LAS_plot(SurveysSample_full)
            ioclass.plot_cummulative_flux(SurveysSample_full)
            exit()

        if 1 == 2:

            surlist = []
            for surveyname in surveynames:
                print('%s/%s' % (folder, surveyname))
                try:
                    survey = iom.unpickleObject('%s/%s' % (folder, surveyname))
                    survey.set_seed_dropout(None)
                    print(survey.relic_filter_kwargs)
                    #survey.relic_filter_kwargs = {"Filter": True, 'minrms': 8, "shape_pca": False}
                    survey.cluster_filter_kwargs = {'minrel': 1, 'zmin': 0.05}
                    for gcl in survey.GCls:
                        for relic in gcl.relics:
                            relic.corrupt_alpha()
                            relic.fpre.label = '$f_\mathrm{Pre}$'

                    surlist.append(survey)
                except:
                    continue

            newdist = lambda x: dbc.measurand(x.Dproj_pix() / x.GCl.R200(), 'Dproj', label='$D_\mathrm{proj,rel}$',
                                              un='$R_{200}$')
            ioclass.create_scattermatrix(surlist,
                                         [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist,
                                          lambda x: x.fpre],
                                          logs=[True, True, False, False, True], suffix='_testPhD', shade=True,
                                          statistics=True) #lambda x: x.GCl.z, lambda x: x.GCl.M200, lambda x: x.LAS,

            exit()

        if 1 == 2:
            distances_list = []
            for surveyname in surveynames:
                print('%s/%s' % (folder, surveyname))
                try:
                    survey = iom.unpickleObject('%s/%s' % (folder, surveyname))
                except:
                    continue

                # just a tweek for now
                survey.relic_filter_kwargs.update({"Filter": True, 'minrms': 8, "shape_pca": True})
                #try:
                #    surmodel = survey.surmodel
                #except:
                #    surmodel = sur.SurModel()
                print(survey.relic_filter_kwargs)
                #survey.set_surmodel(surmodel)
                survey.set_seed_dropout()
                survey.dinfo = survey.GCls[0].dinfo
                survey.scatterkwargs = {"alpha": 0.15, "fmt": "^", "markersize": 7}
                survey.cnt_levels = [9e-4, 1.8e-3, 3.6e-3, 7.2e-3, 1.44e-2, 2.88e-2, 5.76e-2]
                survey.histkwargs = {"alpha": 0.15}
                survey.set_binning(Histo)
                survey.emi_max = 2e-2

                import inference.surveymetrics as surveymetrices
                #metrices = ['number', "flux_kolm", "2DKS", "polarHisto", "polarHisto_simple", "PCA", "alpha"]
                metrices = ["number_cluster", 'number', "polarHisto", "alpha", "2DKS", "flux_kolm"]
                distances = surveymetrices.ABC_dist_severalMetrices(NVSSsurvey, survey, metrics=metrices,
                                         outpath='', delal=False, verbose=True, stochdrop=True)
                print("distances", distances)
                distances_list.append(distances)

            data = np.asarray(distances_list)
            print(data.shape)
            import pandas as pd
            dataset = pd.DataFrame({'folder': folder, 'number_cluster': data[:, 0], 'number': data[:, 1],
                                    'polarHisto': data[:, 2], "alpha": data[:, 3], "2DKS": data[:, 4],
                                    "flux_kolm": data[:, 5]})

            from os.path import normpath, basename
            dataset.to_csv('/data/metric_estimate_%s.csv' % outname)
            np.save("/data/metric_estimate%s" % outname, data)
            exit()


        if 1 == 2:
            SurveysSample = []
            for surveyname in surveynames:
                surveypath = '%s/%s' % (folder, surveyname)
                print(surveypath)
                try:
                    survey = iom.unpickleObject(surveypath)
                except:
                    continue
                #print(survey.name,  [np.log10(eff) for eff in survey.Rmodel.effList],  "np.log10(survey.Rmodel.B0):",
                #                     'np.log10(survey.Rmodel.B0): %.2f' % np.log10(survey.Rmodel.B0), "survey.Rmodel.kappa", survey.Rmodel.kappa)
                #exit()
                survey.seed_dropout = None

                survey.dinfo = survey.GCls[0].dinfo
                survey.scatterkwargs = {"alpha": 0.15, "fmt": "^", "markersize": 7}
                survey.cnt_levels = [9e-4, 1.8e-3, 3.6e-3, 7.2e-3, 1.44e-2, 2.88e-2, 5.76e-2]
                survey.relic_filter_kwargs.update({"Filter": True, 'minrms': 8, "shape_pca": True})
                survey.histkwargs = {"alpha": 0.15}
                survey.set_binning(Histo)
                survey.emi_max = 2e-2
                survey.seed_dropout = None
                survey.set_seed_dropout()
                SurveysSample.append(survey)
            iom.pickleObject_old(SurveysSample, "/data/ManySurveys", append = False)
            print("I did it!!!!")

            SurveysSample = [] #iom.unpickleObject("/data/ManySurveys")
            yeaahh = [NVSSsurvey] + SurveysSample
            print(len(yeaahh))
            ioclass.plot_cummulative_flux(yeaahh)
            exit()

        for surveyname in surveynames:
            SurveysSample = []
            surveypath = '%s/%s' % (folder, surveyname)
            print(surveypath)
            try:
                survey = iom.unpickleObject(surveypath)
            except:
                continue
            #print(survey.name,  [np.log10(eff) for eff in survey.Rmodel.effList],  "np.log10(survey.Rmodel.B0):",
            #                     'np.log10(survey.Rmodel.B0): %.2f' % np.log10(survey.Rmodel.B0), "survey.Rmodel.kappa", survey.Rmodel.kappa)
            #exit()
            survey.seed_dropout = None

            survey.dinfo = survey.GCls[0].dinfo
            survey.scatterkwargs = {"alpha": 0.15, "fmt": "^", "markersize": 7}
            survey.cnt_levels = [9e-4, 1.8e-3, 3.6e-3, 7.2e-3, 1.44e-2, 2.88e-2, 5.76e-2]
            try:
                survey.relic_filter_kwargs.update({"Filter": True, 'minrms': 8, "shape_pca": True})
            except:
                survey.relic_filter_kwargs = {"Filter": True, 'minrms': 8, "shape_pca": True}
            survey.histkwargs = {"alpha": 0.15}
            survey.set_binning(Histo)
            survey.emi_max = 2e-2
            survey.seed_dropout = None
            survey.set_seed_dropout()
            SurveysSample.append(survey)

            #clusters = survey.FilterCluster(**survey.cluster_filter_kwargs)
            #relics = survey.fetch_totalRelics()
            #print(len(survey.GCls), len(clusters), len(relics))

            for gcl in survey.GCls:
                for relic in gcl.relics:
                    relic.corrupt_alpha()




            if 1 == 2:
                usethem = survey.GCls  # survey.FilterCluster()
                Nrelics = np.array([len(gcl.filterRelics(minrms=8, shape=False)) for gcl in usethem])
                M200 = np.array([gcl.M200() for gcl in usethem])
                z = np.array([gcl.z() for gcl in usethem])
                Nrelics = np.expand_dims(Nrelics, axis=1)
                M200 = np.expand_dims(M200, axis=1)
                z = np.expand_dims(z, axis=1)
                measures = np.concatenate((M200, z, Nrelics), axis=1)
                np.save("/data/clusters_with_mass_redshift", measures)



                usethem = survey.GCls  # survey.FilterCluster()
                Nrelics = np.array([len(gcl.filterRelics(minrms=8, shape=False, shape_pca=True)) for gcl in usethem])
                M200 = np.array([gcl.M200() for gcl in usethem])
                z = np.array([gcl.z() for gcl in usethem])
                Nrelics = np.expand_dims(Nrelics, axis=1)
                M200 = np.expand_dims(M200, axis=1)
                z = np.expand_dims(z, axis=1)
                measures = np.concatenate((M200, z, Nrelics), axis=1)
                np.save("/data/clusters_with_mass_redshift_filtered_by_relic", measures)

                survey.set_seed_dropout()
                usethem = survey.FilterCluster(**survey.cluster_filter_kwargs)
                Nrelics = np.array([len(gcl.filterRelics(minrms=8, shape=False)) for gcl in usethem])
                M200 = np.array([gcl.M200() for gcl in usethem])
                z = np.array([gcl.z() for gcl in usethem])
                Nrelics = np.expand_dims(Nrelics, axis=1)
                M200 = np.expand_dims(M200, axis=1)
                z = np.expand_dims(z, axis=1)
                measures = np.concatenate((M200, z, Nrelics), axis=1)
                np.save("/data/clusters_with_mass_redshift_filtered_by_cluster", measures)
                exit()
            survey.set_seed_dropout()





            if 1 == 1:
                import clusterbuster.fitsut as fiut
                # Only the brightest objects
                kwargs = {"minimumLAS": 4, "GClflux": 5}
                usethem = survey.FilterCluster(**kwargs)
                print('len(usethem)', len(usethem), len(survey.FilterCluster(**kwargs)),
                      len(survey.FilterCluster()), len(survey.GCls))
                fiut.sparse_array_to_fits(usethem, survey.outfolder)
                fiut.sparse_array_to_fits(usethem, survey.outfolder, maptype="Subtracted", source_type="compacts")
                fiut.sparse_array_to_fits(usethem, survey.outfolder, maptype="Temperature")
                fiut.sparse_array_to_fits(usethem, survey.outfolder, maptype="Density")
                fiut.sparse_array_to_fits(usethem, survey.outfolder, maptype="Mach")


                #        ioclass.plot_Clusters(survey, dynamicscale=False, relicregions=False, DS9regions=False, sectors=True, colorbar=True)
                ioclass.plot_Clusters(survey, xray=True, dynamicscale=False, relicregions=False, DS9regions=False,
                                      sectors=False, colorbar=True, infolabel=False, subtracted=True,
                                      filterargs={'zmin': 0.05, 'minimumLAS': 0, 'GClflux': 3.6,
                                                  'index': None})
            exit()

            """ DEVELOPMENT 
            if ( -0.9 < survey.Rmodel.kappa < -0.2) and (1.6 > np.log10(survey.Rmodel.B0) > 0.8)  and (-5.0 > np.log10(survey.Rmodel.effList[0]) > -5.3):
                print('___ Choosen:', survey.name)
            else:
                continue
             DEVELOPMENT """
            ioclass.plot_RelicEmission_polar(survey, compsurvey=NVSSsurvey, additive=True, single=False, aligned=True,
                                             cbar=False, mirrored=True, title=None, dpi=600, Histo=Histo)

            """ NVSS """
            SurveysSample.append(NVSSsurvey)

            ioclass.create_shape_LAS_plot([NVSSsurvey, survey])
            ioclass.plot_cummulative_flux([NVSSsurvey, survey])

            ioclass.create_scattermatrix([NVSSsurvey, survey],
                                         [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, lambda x: x.Dproj_pix,
                                          lambda x: x.GCl.z, lambda x: x.GCl.M200, lambda x: x.LAS],
                                            logs=[True, True, False, True, False, True, True], suffix='_test', shade=False)

            """ Replaced density with upstream density """
            ioclass.create_scattermatrix([survey], [lambda x: x.Rho_up, lambda x: x.Mach, lambda x: x.T,
                                                    lambda x: x.Area, lambda x: x.P_rest], suffix='_physical_params')
            #ioclass.create_scattermatrix([survey], [lambda x: x.Rho_up, lambda x: x.Mach, lambda x: x.T,
            #                                        lambda x: x.Area, lambda x: x.P_rest, lambda x: x.fpre],
            #                             suffix='_physical_params_preexisting')

            newdist = lambda x: dbc.measurand(x.Dproj_pix() / x.GCl.R200(), 'Dproj', label='$D_\mathrm{proj,rel}$',
                                              un='$R_{200}$')
            logs_alpha = [True, True, False, True, False, True]
            logs_alpha_new = [True, True, False, False]
            ioclass.create_scattermatrix([NVSSsurvey, survey], [lambda x: x.alpha, lambda x: x.Dproj_pix],
                                         logs=[False, True], shade=False)
            ioclass.create_scattermatrix([NVSSsurvey, survey],
                                         [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist],
                                         logs=logs_alpha, suffix='_newdist', shade=False)
            ioclass.create_scattermatrix([NVSSsurvey, survey],
                                         [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist],
                                         logs=logs_alpha_new, suffix='_newdistB', shade=False)
            ioclass.create_scattermatrix([NVSSsurvey, survey],
                                         [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist,
                                          lambda x: x.GCl.z, lambda x: x.GCl.M200], logs=logs_alpha, suffix='_OWH',
                                          shade=False)
            ioclass.create_scattermatrix([NVSSsurvey, survey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha,
                                                                lambda x: x.Dproj_pix, lambda x: x.GCl.z,
                                                                lambda x: x.GCl.M200], logs=logs_alpha, suffix='_OWH2',
                                         shade=False)
            ioclass.create_scattermatrix([NVSSsurvey, survey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.GCl.z],
                                         logs=logs_alpha, suffix='_z', shade=False)



#                ioclass.create_Power_LLS(     SurveysSample, (zrange, colors, zsymbols), minrel=1)
     
        print('#Surveys in SurveysSample', len(SurveysSample))


    """ Cosmic Web """
    if create_CWpaperplots:
#        ioclass.create_scattermatrix( [NVSSsurvey,NVSSsurvey], (zrange, colors, zsymbols)) 
        ii = 0
        import glob
        for surveypath in glob.glob('%s/CW_02_*' % (folder)):  #['CW_02_5.770','CW_02_5.270','CW_02_5.020','CW_02_4.770','CW_02_4.520','CW_02_4.270','CW_02_4.520']: #'CW_02', 'CW_02_flat','CW_02_modelB_flat','CW_02_modelB']:        
            print('A')
            SurveysSample = []
#            survey = iom.unpickleObject('%s/%s/pickled/Survey' % (folder, surveyname))
            try:
                survey = iom.unpickleObject('%s/pickled/Survey' % (surveypath))
            except:
                continue
            print('B')
            survey.dinfo = survey.GCls[0].dinfo       
            survey.scatterkwargs = {"alpha": 0.15, "fmt": "^", "markersize": 7}
            survey.histkwargs = {"alpha": 0.15}
            survey.hist_main = Histo
            SurveysSample.append(survey)

            
            ioclass.plot_RelicEmission_polar(survey, compsurvey=NVSSsurvey, additive=True, single=False, aligned=True, cbar=False, mirrored=True, title=None, dpi=600, Histo=Histo)
              
#            SurveysSample.append(NVSSsurvey) 
            if len(survey.FilterCluster(**survey.cluster_filter_kwargs)) > 0:
                """ Replace density with upstream density """
                ioclass.create_scattermatrix([survey], [lambda x: x.Rho, lambda x: x.Mach, lambda x: x.T,
                                                        lambda x: x.Area, lambda x: x.P_rest], suffix='_physical_params')

                newdist = lambda x: dbc.measurand(x.Dproj_pix() / x.GCl.R200(), 'Dproj', label='$D_\mathrm{proj,rel}$',
                                                  un='$R_{200}$')
                newalpha = lambda x: dbc.measurand(np.abs(x.alpha()), '|alpha|', label='$|\\alpha|$', un=None)
                logs_alpha = [True, True, False, True, False, True]
                logs_alpha_new = [True, True, False, False]
                ioclass.create_scattermatrix([NVSSsurvey, survey], [lambda x: x.alpha, lambda x: x.Dproj_pix],
                                             logs=[False, True])
                ioclass.create_scattermatrix([NVSSsurvey, survey],
                                             [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist],
                                             logs=logs_alpha, suffix='_newdist')
                ioclass.create_scattermatrix([NVSSsurvey, survey],
                                             [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist],
                                             logs=logs_alpha_new, suffix='_newdistB')
                ioclass.create_scattermatrix([NVSSsurvey, survey],
                                             [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha, newdist,
                                              lambda x: x.GCl.z, lambda x: x.GCl.M200], logs=logs_alpha, suffix='_OWH')
                ioclass.create_scattermatrix([NVSSsurvey, survey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.alpha,
                                                            lambda x: x.Dproj_pix, lambda x: x.GCl.z,
                                                            lambda x: x.GCl.M200], logs=logs_alpha, suffix='_OWH2')
                ioclass.create_scattermatrix([NVSSsurvey, survey], [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.GCl.z],
                                             logs=logs_alpha, suffix='_z')

    return 0



if __name__ == "__main__":
    main()