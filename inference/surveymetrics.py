#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:15:19 2017

@author: jakobg
"""
from __future__ import division, print_function

import clusterbuster.dbclasses      as dbc
import clusterbuster.surveyclasses  as cbclass
#import clusterbuster.surveyut       as suut
 
import os
import glob
import numpy  as np
import ndtest as KSmaster
import math
import time
      
SURVEYCOUNT = 0   # This is the global count for the ABC process initiated, somehow hacky

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
    



def ABC_dist_music2( SurveyA, SurveyB, delal=False, stochdrop=True):
    ''' 
    Returns the distance within the MUSIC-2/NVSS metric
    '''

    for eff in SurveyB.Rmodel.effList:
        
        A,B,C  = ABC_dist_severalMetrices( SurveyA, SurveyB, delal=delal, stochdrop=stochdrop)
        sum_deviation = A+B+C
        print('surveymetrics::ABC_dist_music2()::', sum_deviation, eff, 'A %.2e  B %.2e C %.2e' % (A,B,C))
        returnarg = sum_deviation
    
    ''' We assume that only one model is tested '''
    return returnarg



def ABC_dist_severalMetrices( SurveyA, SurveyB, delal=True, verbose=False, stochdrop=True):
    ''' 
    Returns the distance within the MUSIC-2/NVSS metric
    you have: data,model
    
    SurveyA: realworld
    SurveyB: model
    
    '''
#    
#    ''' DEBUGGIN G:BEGIN'''
#    surveypath = '/data/ClusterBuster-Output/%s' % (SurveyB)
#    with open("%s.txt" % (surveypath), "a") as f:
#        s = np.random.normal(5, 0.3, 3)
#        f.write("Test %.3f %.3f %.3f \n" %(s[0], s[1], s[2]))
#    return s
#    '''DEBUGGING:END'''

    if verbose: print(SurveyA.name, SurveyB.name)
    

    for eff in SurveyB.Rmodel.effList:

        if stochdrop: SurveyB.set_dropseed()
        if len([gcl.updateInformation() for gcl in SurveyB.FilterCluster(minrel=1)]) < 3:
            A,B,C,D   = 1e2,1e2,1e2,1e9

        else:
            SurveyA.FilterCluster(minrel=1)
            SurveyB.FilterCluster(minrel=1)
            print ('SurveyB.filteredClusters', len(SurveyB.GCls), len(SurveyB.filteredClusters))
            A  = ABC_summaryStatistics_numbers([SurveyA,SurveyB])
            B  = ABC_summaryStatistics_polarHisto([SurveyA,SurveyB], eff)
#            C  = 1.0 #ABC_summaryStatistics_2DKS([SurveyA,SurveyB], eff, parA=lambda x: x.largestLAS, parB = lambda y: y.P_rest)
            D  = ABC_summaryStatistics_logMach([SurveyA,SurveyB])
            
            
            ''' Heavyly inspired by https://plot.ly/ipython-notebooks/principal-component-analysis/ '''
            newdist      =  lambda x:  dbc.measurand( x.Dproj_pix()/x.GCl.R200(), 'Dproj',label='$D_\mathrm{proj,rel}$',  un = '$R_{200}$' )
            plotmeasures = [lambda x: x.LLS, lambda x: x.P_rest, lambda x: x.Mach, newdist]
                
                
            from sklearn.preprocessing import StandardScaler
            X = SurveyB.fetchpandas(plotmeasures, surname=False).dropna().as_matrix()   #.data()   #, kwargs_FilterCluster={}, kwargs_FilterObjects={}
            X_std = StandardScaler().fit_transform(X)
                
            
            from sklearn.decomposition import PCA as sklearnPCA
            
            
            '''
            # Make a list of (eigenvalue, eigenvector) tuples
            cor_mat1 = np.corrcoef(X_std.T)
            
            eig_vals, eig_vecs = np.linalg.eig(cor_mat1)
            
            print('Eigenvectors \n%s' %eig_vecs)
            print('\nEigenvalues \n%s' %eig_vals)

            eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]
            
            # Sort the (eigenvalue, eigenvector) tuples from high to low
            eig_pairs.sort()
            eig_pairs.reverse()
            
            # Visually confirm that the list is correctly sorted by decreasing eigenvalues
            print('Eigenvalues in descending order:')
            for i in eig_pairs:
                print(i[0])
            
            matrix_w = np.hstack((eig_pairs[0][1].reshape(4,1), 
                                  eig_pairs[1][1].reshape(4,1)))
            
            print('Matrix W:\n', matrix_w)
            
            Y = X_std.dot(matrix_w)
            '''
            
            sklearn_pca = sklearnPCA(n_components=2)
            Y_sklearn = sklearn_pca.fit_transform(X_std)
            
            ''' This gives you an proxy for the average summed square error in the 2-D dimensional reduction via pca '''
            C = np.sum(Y_sklearn**2)/len(Y_sklearn[0])
            
            
        print('surveymetrics::ABC_dist_severalMetrices::', SurveyA.name, 'VS', SurveyB.name, ' metric disimilarity: A', A, 'B', B, 'C', C, 'D', D)
            
        
        ''' This deletes the survey folder z-snapshot files
        '''
        if delal:
            outpath      = '/data/ClusterBuster-Output/'

            file_path = "%s/pickled/Survey.pickle" % (SurveyB.outfolder)        
            if verbose:
                print('surveymetrics::ABC_dist_severalMetrices:: SurveyB.name, surveyB.outfolder:', SurveyB.name, SurveyB.outfolder)
                print('surveymetrics::ABC_dist_severalMetrices:: file_path, os.path.isfile(file_path)', file_path, os.path.isfile(file_path))
            if os.path.isfile(file_path):             
                n = 0
                while n < 10:
                    try:
                        with open('/data/ClusterBuster-Output/AllSurveys/count.txt', 'r') as f:
                            SURVEYCOUNT = int(f.readline())
                        print(SURVEYCOUNT) # Is implemented to write the correct survey output files for the ABC abbroach
                        os.system("cp -rf %s/pickled/Survey.pickle %s/Survey%05i.pickle" % (SurveyB.outfolder, os.path.join(SurveyB.outfolder, '..', 'AllSurveys'), SURVEYCOUNT))
                        with open('/data/ClusterBuster-Output/AllSurveys/count.txt', 'w') as f:
                            f.write(str(SURVEYCOUNT+1))
                        n=10
                    except:
                        n += 1
                        time.sleep(0.5)
                        print('surveymetrics::ABC_dist_severalMetrices:: Could not write counter.')

                os.system("rm -rf %s/pickled/Survey.pickle" % (SurveyB.outfolder))
                    
                    
            ''' In the future it might be wise to save every Survey output ... you might want to have the walker id and the iteration id,
                I just don't know where to get them from within the run.            
            '''
            
            ''' We increment the current logfile number by one ... just to show how much we have progressed '''
            with open("%s//%s.txt" % (outpath,SurveyB.name), "a") as f:
                Rm  = SurveyB.Rmodel
                eff = SurveyB.Rmodel.effList[0]
                if  A<100:
                    line = "%8.5f %8.5f %8.5f, %8.5f" % (A, B, C, D) 
                else: #DEBUGGING purposses
                    line = "None None None None"
                    
                if isinstance(Rm, cbclass.PreModel_Hoeft):
                    line += ' %+.4e %+.4e %+.4e' % (eff, Rm.B0, Rm.kappa) + ' %+.4e +.4e %+.4e %+.4e %+.4e\n' % (Rm.kappa, Rm.t0, Rm.t1, Rm.n0, Rm.n1)
                if isinstance(Rm, cbclass.PreModel_Gelszinnis):
                    line += ' %+.4e %+.4e %+.4e' % (eff, Rm.B0, Rm.kappa) + ' %+.4e %+.4e %+.4e %+.4e\n' % (Rm.p0, Rm.p_sigma, Rm.sigmoid_0, Rm.sigmoid_width)
                f.write(line)
    

    '''PhD Thesis plot '''
    return [A]     
    '''PhD Thesis END '''
    return [A,B,C,D]
    



def abcpmc_dist_severalMetrices( SurveyA, SurveyB, delal=True, stochdrop=True):
  ''' abcpmc differs from astroABC in that sense that the model data is called before the data in the metric arguement,
      so it is metric(model,data) instead of metric(data,model)
      
  So this is just a wrapper    
  '''
  return ABC_dist_severalMetrices( SurveyB, SurveyA, delal=True, stochdrop=stochdrop)