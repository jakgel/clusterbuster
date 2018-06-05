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
import numpy  as np
import ndtest as KSmaster
import math
import time

'''=============== Baustelle: Imlement in Metric & Run Survey'''
def Clusters_discovery_prop(survey, discovery_prop=None, maxcomp=None, verbose=False):
    ''' Return a weighted relic number count within all clusters
    input: either  a Clusterbuster Class Surveys
           or      a GalaxyClusterLists including relics
           and the efficiency at which to count the relics
           optional: the function of the discovery prop that is a shape&LAS dependent discovery probability
    It allows to filter for cetain shape regimes      
    Returns
    ------
    distance: float
    '''
    
#    import scipy.stats.mstats as mstats
    if isinstance(survey, cbclass.Survey):
        ''' Assume to work with surveys '''
        relics = survey.fetch_totalRelics(maxcomp=maxcomp)
    else:
        ''' Asume to work with lists of galaxy clusters '''
        relics = [gcl.filterRelics(maxcomp=maxcomp)          for gcl in survey]

    weightedsum = len(relics)
    if verbose:
        print('surveymetrics::Clusters_discovery_prop()::Survey:', survey.name, weightedsum)

    return weightedsum
'''==============='''


def ABC_summaryStatistics_polarHisto(Surveys, eff):
    ''' Compares the 'average relic' of survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
           or      2 Polary binned Histogramms of the same dimensions
           
           
    This method will fail if the first survey class doesn't have any relics!
    '''
    [SurveyA,SurveyB] = Surveys

    # Filtering by redshift
    zborder = 0.05
    ztype   = '>'
    
    norm = dbc.norm('R200',Nexp=1.5) # I used Nexp=2.0 in the past! ANd maybe I just should use 1.0!!!!!!!!!!!!!!
    for Survey in Surveys:
        Survey.mainhist = dbc.Histogram2D(nbins=(32,30), fromto= [[0,2.*np.pi],[0,1.5]], norm=norm )     # angle_projected(rad), D_proj(R200) 
        Survey.expScale = 0.65
    
    if SurveyB.polar()[0] is None :
        HiB = 0
    else:
        HiB = SurveyB.polar(zborder = zborder,ztype = ztype, normalize=True)[1][0]
    HiA = SurveyA.polar(zborder = zborder,ztype = ztype, normalize=True)[1][0]
    deviation =  np.sum(np.abs(HiA-HiB))
    
    return deviation
    

def ABC_summaryStatistics_2DKS(Surveys, eff, parA=lambda x: x.M200, parB = lambda y: y.P_rest):
    ''' Compares survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
    '''
    
    [A, B] = Surveys
    
#    import scipy.stats.mstats as mstats
    if isinstance(A, cbclass.Survey) and isinstance(B, cbclass.Survey):
        ''' Assume to work with surveys '''
        cl_A = [gcl.updateInformation()        for gcl in A.GCls]
        cl_B = [gcl.updateInformation(eff=eff) for gcl in B.filteredClusters]
    else:
        ''' Asume to work with lists of galaxy clusters '''
        cl_A = [gcl.updateInformation()        for gcl in A]
        cl_A = [gcl.updateInformation(eff=eff) for gcl in B]
        
    x1  = np.asarray([parA(gcl).value for gcl in cl_A])
    y1  = np.asarray([parB(gcl).value for gcl in cl_A])   

    x2  = np.asarray([parA(gcl).value for gcl in cl_B])
    y2  = np.asarray([parB(gcl).value for gcl in cl_B])   
    
    
    '''Two-dimensional Kolmogorov-Smirnov test on two samples. 
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
    '''

    if x1.shape[0] > 2 and x2.shape[0] > 2:
        access, stats = KSmaster.ndtest.ks2d2s(x1, y1, x2, y2, nboot=None, extra=True)
        if math.isnan(access): access = 1.e-10
    else:
        access = 1.e-10
    print('surveymetrics::ABC_summaryStatistics_2DKS()::access', access)
    return np.log10(1/access)

def ABC_summaryStatistics_numbers(Surveys, maxcomp=None, verbose = False):
    ''' Compares survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
           or      2 GalaxyClusterLists including relics
           and the efficiency at which to compare the galaxy clusters
           optional: the function of the discovery prop that is a shape dependent discovery probability
    It allows to filter for cetain shape regimes      
    Returns:
        All radio relics
    ------
    distance: float
    '''
    [A, B] = Surveys

    sum_A = Clusters_discovery_prop(A, maxcomp=maxcomp)
    sum_B = Clusters_discovery_prop(B, maxcomp=maxcomp)

    if verbose:
        print('surveymetrics::ABC_summaryStatistics_numbers::Survey:', A.name, len(A.GCls)            , sum_A)
        print('surveymetrics::ABC_summaryStatistics_numbers::Model :', B.name, len(B.filteredClusters), sum_B)
    
    # This is like assuming a students distribution (?) in the number count and taking the the reduced sqrt of the number as the standart deviation
    #deviation =  np.abs(sum_A-sum_B)  / max(1.,np.sqrt(sum_A - 1.))
    deviation =  np.abs(sum_A-sum_B)  / np.sqrt(sum_A*max(1.,sum_B))
    return deviation



def ABC_summaryStatistics_logMach(Surveys):
    ''' Compares survey B (simulation) with survey A (real world survey)
        Derives the Histogram of the log(Mach) number and returns the differnce of the medians
    '''
    [A, B] = Surveys

    if isinstance(A, cbclass.Survey):
        ''' Assume to work with surveys '''
        relicsA = A.fetch_totalRelics()
        relicsB = B.fetch_totalRelics()
    else:
        ''' Asume to work with lists of galaxy clusters '''
        relicsA = [gcl.filterRelics() for gcl in A]
        relicsB = [gcl.filterRelics() for gcl in B]      
    
    
    # Get mach-numbers and remove nans
    A = np.array( [np.log10(relic.Mach()) for relic in relicsA]) 
    B = np.array( [np.log10(relic.Mach()) for relic in relicsB]) 
    
    mA = A[~np.isnan(A)].mean()
    mB = B[~np.isnan(B)].mean()
    
    return abs(mA-mB)
    

def ABC_dist_severalMetrices( SurveyA, SurveyB,  metrics = ['numbers'],
                             outpath = '', delal=True, verbose=False, stochdrop=True):
    ''' 
    Returns the distance within the MUSIC-2/NVSS metric
    you have: data,model
    
    SurveyA: realworld
    SurveyB: model
    
    '''
    print('ABC_dist_severalMetrices',metrics)

    if verbose: print(SurveyA.name, SurveyB.name)
    

    ''' The efficiency is outdated and should be replaced asap '''
    for eff in SurveyB.Rmodel.effList:

        if stochdrop: SurveyB.set_dropseed()
        
        if len([gcl.updateInformation() for gcl in SurveyB.FilterCluster(minrel=1)]) < 3:
            distances = [1e9 for m in metrics]
        else:  
            distances = []
            SurveyA.FilterCluster(minrel=1)
            SurveyB.FilterCluster(minrel=1)
            print ('SurveyB.filteredClusters', len(SurveyB.GCls), len(SurveyB.filteredClusters))
            for metric in metrics:
                print('metric', metric)
                
                if metric == 'number':
                    distance  = ABC_summaryStatistics_numbers([SurveyA,SurveyB])
                if metric == 'polarHisto':
                    distance  = ABC_summaryStatistics_polarHisto([SurveyA,SurveyB], eff)
                if metric == 'logMach':
                    distance  = ABC_summaryStatistics_logMach([SurveyA,SurveyB])
                if metric == 'PCA':       
                    ''' Heavily inspired by https://plot.ly/ipython-notebooks/principal-component-analysis/ '''
                    newdist      =  lambda x:  dbc.measurand( x.Dproj_pix()/x.GCl.R200(), 'Dproj',label='$D_\mathrm{proj,rel}$',  un = '$R_{200}$' )
                    plotmeasures = [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.Mach, newdist]
                    
                    
                    from sklearn.preprocessing import StandardScaler
                    X = SurveyB.fetchpandas(plotmeasures, surname=False).dropna().as_matrix()   #.data()   #, kwargs_FilterCluster={}, kwargs_FilterObjects={}
                    X_std = StandardScaler().fit_transform(X)
                        
                    
                    from sklearn.decomposition import PCA as sklearnPCA
                    
                    
                    sklearn_pca = sklearnPCA(n_components=2)
                    Y_sklearn = sklearn_pca.fit_transform(X_std)
                    
                    ''' This gives you an proxy for the average summed square error in the 2-D dimensional reduction via pca '''
                    distance = np.sum(Y_sklearn**2)/len(Y_sklearn[0])
                
                distances.append(distance)
                     
        print('surveymetrics::ABC_dist_severalMetrices::', SurveyA.name, 'VS', SurveyB.name, ' metric disimilarity:', 'metric disimilarity:', ['%s: %.3e' % (m,d) for (m,d) in zip(metrics,distances)] )
            
        
        ''' This puts the survey to a bagged sutvey folde rand increases the counter. It might be interesting to know if this number is also the number of the runs.
        '''
        if delal:

            file_path = "%s/pickled/Survey.pickle" % (SurveyB.outfolder)        
            if verbose:
                print('surveymetrics::ABC_dist_severalMetrices:: SurveyB.name, surveyB.outfolder:', SurveyB.name, SurveyB.outfolder)
                print('surveymetrics::ABC_dist_severalMetrices:: file_path, os.path.isfile(file_path)', file_path, os.path.isfile(file_path))
            if os.path.isfile(file_path):             
                n = 0
                while n < 10:  
                    try:
                        with open('%s/count.txt' % (outpath), 'r') as f:
                            SURVEYCOUNT = int(f.readline())
                        print(SURVEYCOUNT) # Is implemented to write the correct survey output files for the ABC abbroach
                        with open('%s/count.txt' % (outpath), 'w') as f:
                            f.write(str(SURVEYCOUNT+1))
                        iom.check_mkdir(outpath+'surveys/') 
                        os.system("cp -rf %s/pickled/Survey.pickle %s/Survey%05i.pickle" % (SurveyB.outfolder, os.path.join(SurveyB.outfolder, '..', 'AllSurveys'), SURVEYCOUNT))  
                        n=10
                    except:
                        n += 1
                        time.sleep(0.5)
                        print('surveymetrics::ABC_dist_severalMetrices:: Could not write counter.')
                        print("___ cp -rf %s/pickled/Survey.pickle %s/Survey_unknown.pickle" % (SurveyB.outfolder, os.path.join(SurveyB.outfolder, '..', 'AllSurveys')))

                os.system("rm -rf %s/pickled/Survey.pickle" % (SurveyB.outfolder))
                    
                
            
            ''' We increment the current logfile number by one ... just to show how much we have progressed '''
            with open("%s/logfile.txt" % (outpath), "a") as f:
                Rm  = SurveyB.Rmodel
                eff = SurveyB.Rmodel.effList[0]

                line = ''
                for dist in distances:
                    line += "%8.5e " % (dist)
                line += '%+.4e %+.4e %+.4e' % (eff, Rm.B0, Rm.kappa)
      
                if isinstance(Rm, cbclass.PreModel_Hoeft):
                    line += ' %+.4e +.4e %+.4e %+.4e %+.4e\n' % (Rm.kappa, Rm.t0, Rm.t1, Rm.n0, Rm.n1)
                if isinstance(Rm, cbclass.PreModel_Gelszinnis):
                    line += ' %+.4e %+.4e %+.4e %+.4e\n' % (Rm.p0, Rm.p_sigma, Rm.sigmoid_0, Rm.sigmoid_width)
                f.write(line)
    
    return distances
    



def abcpmc_dist_severalMetrices( SurveyA, SurveyB, metrics=['number'], outpath = '', delal=True, stochdrop=True):
  ''' abcpmc differs from astroABC in that sense that the model data is called before the data in the metric arguement,
      so it is metric(model,data) instead of metric(data,model)
      
  So this is just a wrapper    
  '''
  return ABC_dist_severalMetrices( SurveyB, SurveyA, metrics=metrics, outpath = outpath, delal=delal, stochdrop=stochdrop)