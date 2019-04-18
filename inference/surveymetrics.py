#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:15:19 2017

@author: jakobg
"""
from __future__ import division, print_function

import clusterbuster.dbclasses      as dbc
import clusterbuster.surveyclasses  as cbclass
import clusterbuster.iout.misc      as iom
#import clusterbuster.surveyut       as suut
 
import os
import shutil
import numpy  as np
import ndtest
import math
import time

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
import scipy

def abcpmc_dist_severalMetrices(SurveyA, SurveyB, metrics=['number'], outpath='', delal=True, stochdrop=True):
  """ abcpmc differs from astroABC in that sense that the model data is called before the data in the metric arguement,
      so it is metric(model,data) instead of metric(data,model)
      
  So this is just a wrapper    
  """
  return ABC_dist_severalMetrices(SurveyB, SurveyA, metrics=metrics, outpath = outpath, delal=delal, stochdrop=stochdrop)


def ABC_dist_severalMetrices(SurveyA, SurveyB, metrics=['number'], outpath='', delal=True,
                             verbose=False, stochdrop=True):
    """ 
    Returns the distance within the MUSIC-2/NVSS metric
    you have: data,model
    
    SurveyA: realworld
    SurveyB: model
    
    """

    print('ABC_dist_severalMetrices', metrics)
    if verbose:
        print(SurveyA.name, SurveyB.name)

    if stochdrop:
        SurveyB.set_seed_dropout()

    distances = []
    SurveyA.FilterCluster(**SurveyA.cluster_filter_kwargs)
    SurveyB.FilterCluster(**SurveyB.cluster_filter_kwargs)

    if len(SurveyB.filteredClusters) < 3:
        distances = [1e9 for m in metrics]
        return distances

    print('SurveyA.GCls', len(SurveyA.GCls), '-->', 'SurveyB.filteredClusters', len(SurveyA.filteredClusters))
    print('SurveyB.GCls', len(SurveyB.GCls), '-->', 'SurveyB.filteredClusters', len(SurveyB.filteredClusters))

    relicsA = SurveyA.fetch_totalRelics()
    A = np.array([min(-1, relic.alpha()) for relic in relicsA])
    for metric in metrics:
        print('metric:', metric)
        if metric == 'number':
            distance = ABC_summaryStatistics_number_relics([SurveyA,SurveyB], verbose=verbose)
        elif metric == 'number_cluster':
            distance = ABC_summaryStatistics_number_cluster([SurveyA,SurveyB], verbose=verbose)
        elif metric == 'flux_kolm':
            distance = ABC_summaryStatistics_flux_komogorov([SurveyA,SurveyB])
        elif metric == 'polarHisto':
            distance = ABC_summaryStatistics_polarHisto([SurveyA,SurveyB])
        elif metric == 'polarHisto_simple':
            distance = ABC_summaryStatistics_polarHisto_simple([SurveyA,SurveyB])
        elif metric == 'logMach':
            distance = ABC_summaryStatistics_logMach([SurveyA,SurveyB])
        elif metric == 'alpha':
            distance = ABC_summaryStatistics_alpha([SurveyA,SurveyB])
        elif metric == '2DKS':
            distance = ABC_summaryStatistics_2DKS([SurveyA,SurveyB])
        elif metric == 'PCA':
            distance = ABC_summaryStatistics_PCA([SurveyA,SurveyB])
        else:
            print(metric, 'is not an implemented metric!')
        distances.append(distance)

    print('surveymetrics::ABC_dist_severalMetrices::', SurveyA.name, 'VS', SurveyB.name, 'metric disimilarity:',
          ['%s: %.3e' % (m,d) for (m,d) in zip(metrics,distances)])


    """ This puts the survey to a bagged survey folder and increases the counter.
    It might be interesting to know if this number is also the number of the runs.
    """
    if delal:
        file_path = "%s/pickled/Survey.pickle" % (SurveyB.outfolder)
        if os.path.isfile(file_path):
            n = 0
            while n < 10:
                try:
                    with open('%s/count.txt' % (outpath), 'r') as f:
                        SURVEYCOUNT = int(float(f.readline()))
                    print('SURVEYCOUNT', SURVEYCOUNT) # Is implemented to write the correct survey output files for the ABC abbroach
                    with open('%s/count.txt' % (outpath), 'w') as f:
                        f.write(str(SURVEYCOUNT+1))

                    SurveyB.name = "%s_%05i" % (SurveyB.name_short, SURVEYCOUNT)
                    outfolder_old = SurveyB.outfolder
                    SurveyB.outfolder = '/data/ClusterBuster-Output/%s' % (SurveyB.name)

                    if verbose:
                        print('surveymetrics::ABC_dist_severalMetrices:: SurveyB.name, surveyB.outfolder:',
                              SurveyB.name, SurveyB.outfolder)
                        print('surveymetrics::ABC_dist_severalMetrices:: file_path, os.path.isfile(file_path)',
                              file_path, os.path.isfile(file_path))

                    iom.check_mkdir(outpath+'surveys/')
                    iom.pickleObject(SurveyB, "%s/surveys/" % outpath, "Survey_%05i" % SURVEYCOUNT)  #obj, location, oname, append = False
                    print("shutil.rmtree:", outfolder_old)
                    shutil.rmtree(outfolder_old)
                    n = 10
                except:
                    n += 1
                    time.sleep(0.2)
                    print('surveymetrics::ABC_dist_severalMetrices:: Could not write counter.')
                    print("___ cp -rf %s/pickled/Survey.pickle %s/surveys/Survey_unknown.pickle" % (SurveyB.outfolder, outpath))

        """ We increment the current logfile number by one ... just to show how much we have progressed """
        with open("%s/logfile.txt" % outpath, "a") as f:
            Rm  = SurveyB.Rmodel
            Sm  = SurveyB.surmodel
            eff = SurveyB.Rmodel.effList[0]

            line = ''
            for dist in distances:
                line += "%8.5e " % dist
            line += '%7i %+.4e %+.4e %+.4e' % (SURVEYCOUNT, eff, Rm.B0, Rm.kappa)

            if isinstance(Rm, cbclass.PreModel_Hoeft):
                line += ' %+.4e %+.4e %+.4e %+.4e %+.4e' % (Rm.t0, Rm.t1, Rm.ratio, Rm.n0, Rm.n1)
            if isinstance(Rm, cbclass.PreModel_Gelszinnis):
                line += ' %+.4e %+.4e %+.4e %+.4e' % (Rm.p0, Rm.p_sigma, Rm.sigmoid_0, Rm.sigmoid_width)
            if Sm is not None:
                line += ' %+.4e %+.4e' % (Sm.relic_filter_pca_a, Sm.relic_filter_pca_b)
            line += '\n'

            f.write(line)
    
    return distances

def ABC_summaryStatistics_number_relics(Surveys, verbose=False):
    """ Compares survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
           or      2 GalaxyClusterLists including relics
           and the efficiency at which to compare the galaxy clusters
           optional: the function of the discovery prop that is a shape dependent discovery probability
    It allows to filter for cetain shape regimes
    Returns:
        All radio relics
    ------
    distance: float
    """
    [A, B] = Surveys

    sum_A = len(A.fetch_totalRelics())
    sum_B = len(B.fetch_totalRelics())

    if verbose:
        print('surveymetrics::ABC_summaryStatistics_number_relics::Survey:', A.name, sum_A)
        print('surveymetrics::ABC_summaryStatistics_number_relics::Model :', B.name, sum_B)

    # This is like assuming a students distribution (?) in the number count and taking the the reduced sqrt of the number as the standart deviation
    # deviation =  np.abs(sum_A-sum_B)  / max(1.,np.sqrt(sum_A - 1.))
    print(sum_A, sum_B, max(1., sum_A), max(1., sum_B))
    deviation = np.abs(sum_A - sum_B) / np.sqrt(max(1., sum_A) * max(1., sum_B))
    return deviation


def ABC_summaryStatistics_number_cluster(Surveys, verbose=False):
    """ Compares survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
           or      2 GalaxyClusterLists including relics
           and the efficiency at which to compare the galaxy clusters
           optional: the function of the discovery prop that is a shape dependent discovery probability
    It allows to filter for cetain shape regimes
    Returns:
        All radio relics
    ------
    distance: float
    """
    [A, B] = Surveys

    sum_A = len(A.filteredClusters)
    sum_B = len(B.filteredClusters)
    if verbose:
        print("A.name", A.name, 'len(A.GCls)', len(A.GCls), 'len(A.filteredClusters)', len(A.filteredClusters))
        print("B.name", B.name, 'len(B.GCls)', len(B.GCls), 'len(B.filteredClusters)', len(B.filteredClusters))
    # This is like assuming a students distribution (?) in the number count and taking the the reduced sqrt of the number as the standart deviation
    # deviation =  np.abs(sum_A-sum_B)  / max(1.,np.sqrt(sum_A - 1.))
    deviation = np.abs(sum_A - sum_B) / np.sqrt(max(1., sum_A) * max(1., sum_B))
    return deviation


"""==============="""


def ABC_summaryStatistics_polarHisto(Surveys):
    """ Compares the 'average relic' of survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
           or      2 Polary binned Histogramms of the same dimensions
           
           
    This method will fail if the first survey class doesn't have any relics!
    """
    [SurveyA, SurveyB] = Surveys

    norm = dbc.norm('R200', Nexp=1.0)         # In the past I used 1.5 ... because this could be a better scaling
    for Survey in Surveys:
        Survey.hist_main = dbc.Histogram2D(nbins=(32,30), fromto=[[0,2.*np.pi], [0,1.5]], norm=norm)     # angle_projected(rad), D_proj(R200)
        Survey.set_binning()
        Survey.expScale = 0.75

    polar = SurveyB.polar(normalize=True)
    if polar[0] is None:
        HiB = 0
    else:
        HiB = polar[1][0]
    HiA = SurveyA.polar(normalize=True)[1][0]
    deviation = np.sum(np.abs(HiA-HiB))
    
    return deviation


def ABC_summaryStatistics_polarHisto_simple(Surveys):
    """ Compares the 'average relic' of survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
           or      2 Polary binned Histograms of the same dimensions


    This method will fail if the first survey class doesn't have any relics!
    """
    [SurveyA, SurveyB] = Surveys
    norm = dbc.norm('R200', Nexp=1.0)
    for Survey in Surveys:
        Survey.mainhist = dbc.Histogram2D(nbins=(32, 30), fromto=[[0, 2. * np.pi], [0, 1.5]],
                                          norm=norm)  # angle_projected(rad), D_proj(R200)
        Survey.expScale = 0.75

    if SurveyB.polar()[0] is None:
        HiB = 0
    else:
        HiB = SurveyB.polar(normalize=True)[1][0]
    HiA = SurveyA.polar(normalize=True)[1][0]
    deviation = np.sum(np.abs(HiA - HiB))

    return deviation



def ABC_summaryStatistics_flux_komogorov(Surveys, par=lambda x: x.flux):

    [A, B] = Surveys
    fluxes_A = [par(relic).value for relic in A.fetch_totalRelics()]
    fluxes_B = [par(relic).value for relic in B.fetch_totalRelics()]
    statistic, p_value = scipy.stats.ks_2samp(fluxes_A, fluxes_B)

    return 1-p_value

def ABC_summaryStatistics_2DKS(Surveys, verbose=False, parA=lambda x: x.LLS, parB = lambda y: y.P_rest, clusterwide=False):
    """ Compares survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
    """
    
    [A, B] = Surveys
    
#    import scipy.stats.mstats as mstats
    if clusterwide:
        entities_A = [gcl.updateInformation() for gcl in A.GCls]
        entities_B = [gcl.updateInformation() for gcl in B.filteredClusters]
    else:
        entities_A = A.fetch_totalRelics()
        entities_B = B.fetch_totalRelics()


    x1 = np.asarray([parA(ent).value for ent in entities_A])
    y1 = np.asarray([parB(ent).value for ent in entities_A])

    x2 = np.asarray([parA(ent).value for ent in entities_B])
    y2 = np.asarray([parB(ent).value for ent in entities_B])
    
    
    """Two-dimensional Kolmogorov-Smirnov test on two samples. 
    Parameters
    ----------
    x1, y1 : ndarray, shape (n1, )
        Data of sample 1.
    x2, y2 : ndarray, shape (n2, )
        Data of sample 2. Size of two samples can be different.
    extra: bool, optional
        If True, KS statistic is also returned. Default is False.

    Returns
    -------
    p : float
        Two-tailed p-value.
    D : float, optional
        KS statistic. Returned if keyword `extra` is True.
    """

    if x1.shape[0] > 2 and x2.shape[0] > 2:
        p_value, stats = ndtest.ks2d2s(x1, y1, x2, y2, nboot=None, extra=True) #KSmaster.
        distance = 1 - p_value #np.log10(1/p_value)
    else:
        distance = 1
    if math.isnan(p_value):
        distance = 1

    if verbose:
        print('surveymetrics::ABC_summaryStatistics_2DKS()::', distance)
    return distance

def ABC_summaryStatistics_logMach(Surveys):
    """ Compares survey B (simulation) with survey A (real world survey)
        Derives the Histogram of the log(Mach) number and returns the difference of the means
    """
    [A, B] = Surveys

    if isinstance(A, cbclass.Survey):
        """ Assume to work with surveys """
        relicsA = A.fetch_totalRelics()
        relicsB = B.fetch_totalRelics()
    else:
        """ Asume to work with lists of galaxy clusters """
        relicsA = [gcl.filterRelics() for gcl in A]
        relicsB = [gcl.filterRelics() for gcl in B]      
    
    
    # Get mach-numbers and remove nans
    A = np.array([np.log10(relic.Mach()) for relic in relicsA])
    B = np.array([np.log10(relic.Mach()) for relic in relicsB])
    
    mA = A[~np.isnan(A)].mean()
    mB = B[~np.isnan(B)].mean()
    
    return abs(mA-mB)



def ABC_summaryStatistics_alpha(Surveys):
    """ Compares survey B (simulation) with survey A (real world survey)
        Derives the difference of the means in alpha
    """
    [A, B] = Surveys

    #for gcl in A.GCls:
    #    for relic in gcl.relics:
    #        relic.corrupt_alpha()


    if isinstance(A, cbclass.Survey):
        """ Assume to work with surveys """
        relicsA = A.fetch_totalRelics()
        relicsB = B.fetch_totalRelics()
    else:
        """ Asume to work with lists of galaxy clusters """
        relicsA = [gcl.filterRelics() for gcl in A]
        relicsB = [gcl.filterRelics() for gcl in B]

    #corrupt the measurement of the spectral index
    for relic in relicsB:
        relic.corrupt_alpha()

    # Get alpha and remove nans
    A = np.array([min(-1, relic.alpha()) for relic in relicsA])
    B = np.array([min(-1, relic.alpha()) for relic in relicsB])

    mA = A[~np.isnan(A)].mean()
    mB = B[~np.isnan(B)].mean()

    return abs(mA - mB)


def ABC_summaryStatistics_PCA(Surveys):
    """ Heavily inspired by https://plot.ly/ipython-notebooks/principal-component-analysis/ """
    [A, B] = Surveys
    
    newdist = lambda x:  dbc.measurand( x.Dproj_pix()/x.GCl.R200(), 'Dproj', label='$D_\mathrm{proj,rel}$', un='$R_{200}$' )
    plotmeasures = [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.Mach, newdist]

    X1 = A.fetch_pandas(plotmeasures, surname=False).dropna().as_matrix()   #.data()   #, kwargs_FilterCluster={}, kwargs_FilterObjects={}
    X2 = B.fetch_pandas(plotmeasures, surname=False).dropna().as_matrix()   #.data()   #, kwargs_FilterCluster={}, kwargs_FilterObjects={}
    
    X1_std = StandardScaler().fit_transform(X1)
    X2_std = StandardScaler().fit_transform(X2)   
    
#    # http://scikit-learn.org/stable/auto_examples/decomposition/plot_pca_iris.html
#    plt.cla()
#    pca = decomposition.PCA(n_components=3)
#    pca.fit(X)
#    X = pca.transform(X)
     
    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(X1_std)
    
    """ This gives you a proxy for the average summed square error in the 2-D dimensional reduction via pca """
    distance = np.sum(Y_sklearn**2)/len(Y_sklearn[0])
    return distance