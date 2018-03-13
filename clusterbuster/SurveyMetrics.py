#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:15:19 2017

@author: jakobg
"""
from __future__ import division, print_function
import pyutil_my.ClusterBuster_pythonClass as CBclass
import pyutil_my.Custom_DatabaseClasses    as cdb 
import numpy as np
import ndtestMaster                        as KSmaster
import math
#import pyutil_my.IOutil                    as iout
      


'''=============== Baustelle: Imlement in Metric & Run Survey'''
def Clusters_discovery_prop(survey, eff=1, discovery_prop=None, maxcomp=None, verbose=False):
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
    if isinstance(survey, CBclass.Survey):
        ''' Assume to work with surveys '''
        relics = survey.fetch_totalRelics(maxcomp=maxcomp, eff=eff)
    else:
        ''' Asume to work with lists of galaxy clusters '''
        relics = [gcl.filterRelics(maxcomp=maxcomp, eff=eff)          for gcl in survey]

    weightedsum = len(relics)
    if verbose:
        print('Survey:', survey.name, weightedsum)

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
    
    norm = cdb.norm('R200',Nexp=1.5) # I used Nexp=2.0 in the past! ANd maybe I just should use 1.0!!!!!!!!!!!!!!
    for Survey in Surveys:
        Survey.mainhist = cdb.Histogram2D(nbins=(32,30), fromto= [[0,2.*np.pi],[0,1.5]], norm=norm )     # angle_projected(rad), D_proj(R200) 
        Survey.expScale = 0.65
    
    if SurveyB.polar()[0] is None :
        HiA = SurveyA.polar(zborder = zborder, ztype = ztype, normalize=True)[1][0]
        deviation  =  np.sum(np.abs(HiA-0))
    else:
        HiA = SurveyA.polar(zborder = zborder,ztype = ztype, normalize=True)[1][0]
        HiB = SurveyB.polar(zborder = zborder,ztype = ztype, normalize=True)[1][0]
        deviation =  np.sum(np.abs(HiA-HiB))
    
    return deviation
    

def ABC_summaryStatistics_2DKS(Surveys, eff, parA=lambda x: x.M200, parB = lambda y: y.P_rest):
    ''' Compares survey B (simulation) with survey A (real world survey)
    input: either  2 Clusterbuster Class Surveys
    '''
    
    [A, B] = Surveys
    
#    import scipy.stats.mstats as mstats
    if isinstance(A, CBclass.Survey) and isinstance(B, CBclass.Survey):
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
    print('ABC_summaryStatistics_2DKS:', access)
    return np.log10(1/access)

def ABC_summaryStatistics_numbers(Surveys, eff, maxcomp=None, verbose = False):
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
    sum_B = Clusters_discovery_prop(B, maxcomp=maxcomp, eff=eff)

    if verbose:
        print('Survey:', A.name, len(A.GCls)            , sum_A)
        print('Model :', B.name, len(B.filteredClusters), sum_B)
    
    # This is like assuming a students distribution (?) in the number count and taking the the reduced sqrt of the number as the standart deviation
    #deviation =  np.abs(sum_A-sum_B)  / max(1.,np.sqrt(sum_A - 1.))
    deviation =  np.abs(sum_A-sum_B)  / np.sqrt(sum_A*max(1.,sum_B))
    return deviation



def ABC_summaryStatistics_logMach(Surveys, eff):
    ''' Compares survey B (simulation) with survey A (real world survey)
        Derives the Histogram of the log(Mach) number and returns the differnce of the medians
    '''
    [A, B] = Surveys

    
    import scipy.stats as stats
    if isinstance(A, CBclass.Survey):
        ''' Assume to work with surveys '''
        relicsA = A.fetch_totalRelics(eff=eff)
        relicsB = B.fetch_totalRelics(eff=eff)
    else:
        ''' Asume to work with lists of galaxy clusters '''
        relicsA = [gcl.filterRelics(eff=eff) for gcl in A]
        relicsB = [gcl.filterRelics(eff=eff) for gcl in B]      
    
    mA = stats.mean([np.log10(relic.mach()) for relic in relicsA])
    mB = stats.mean([np.log10(relic.mach()) for relic in relicsB])
    
    
    return abs(mA-mB)
    



def ABC_dist_music2( SurveyA, SurveyB, delal=False, stochdrop=True):
    ''' 
    Returns the distance within the MUSIC-2/NVSS metric
    '''

    for eff in SurveyB.Rmodel.effList:
        
        A,B,C  = ABC_dist_severalMetrices( SurveyA, SurveyB, delal=delal, stochdrop=stochdrop)
        sum_deviation = A+B+C
        print('ABC_dist_music2:', sum_deviation, eff, 'A %.2e  B %.2e C %.2e' % (A,B,C))
        returnarg = sum_deviation
    
    ''' We assume that only one model is tested '''
    return returnarg



def ABC_dist_severalMetrices( SurveyA, SurveyB, delal=True, verbose=False, stochdrop=True):
    ''' 
    Returns the distance within the MUSIC-2/NVSS metric
    you have: data,model
    
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
        if len([gcl.updateInformation() for gcl in SurveyB.FilterCluster(minrel=1)]) < 1:
            A,B,C   = 1e2,1e2,1e2
#            print([len(gcl.relics)  for gcl in SurveyB.GCls])
#            print('!!!!')
        else:
            print ('SurveyB.filteredClusters', len(SurveyB.GCls), len(SurveyB.filteredClusters))
#            print('A')
            A  = ABC_summaryStatistics_polarHisto([SurveyA,SurveyB], eff)
            B  = ABC_summaryStatistics_numbers([SurveyA,SurveyB], eff)
            C  = ABC_summaryStatistics_2DKS([SurveyA,SurveyB], eff, parA=lambda x: x.largestLAS, parB = lambda y: y.P_rest)
            D  = ABC_summaryStatistics_logMach([SurveyA,SurveyB], eff)
        print(SurveyA.name, SurveyB.name, ': A', A, 'B', B, 'C', C)
            
        
        
        ''' This deletes the survey folder z-snapshot files
        '''
        if delal:
            outpath      = '/data/ClusterBuster-Output/'
            surveypath   = outpath + '%s/' % (SurveyB.name)

#            import shutil    
#            shutil.rmtree(surveypath)     # This is to brutal, we would like to keep at least the survey id
            
            #####
            import os
            import glob
            
            file_path = "%s/pickled/SurveySample.pickle" % (surveypath)        
            if verbose:
                print(surveypath)
                print(file_path)
            if os.path.isfile(file_path):
                #Deletes unneeded files.
                for CleanUp in glob.glob("%s/pickled/*.pickle" % (surveypath) ):
                    if not CleanUp.endswith('SurveySample.pickle'):    
                        os.remove(CleanUp)
                        
            ''' In the future it might be wise to save every Survey output ... you might want to have the walker id and the iteration id,
                I just don't know where to get them from within the run.            
            '''
            
            ''' We increment the current logfile number by one ... just to show how much we have progressed '''
            with open("%s//%s.txt" % (outpath,SurveyB.name), "a") as f:
                Rm  = SurveyB.Rmodel
                eff = SurveyB.Rmodel.effList[0]
                if  A<100:
                    line = "%8.5f %8.5f %8.5f" % (A, B, C) 
                else: #DEBUGGING purposses
                    line = "None None None"
                line += ' %+.4e %+.4e %+.4e' % (eff, Rm.B0, Rm.kappa) + ' %+.4e %+.4e %+.4e %+.4e\n' % (Rm.p0, Rm.p_sigma, Rm.sigmoid_0, Rm.sigmoid_width)
                f.write(line)
            

    return [A,B,C]



def abcpmc_dist_severalMetrices( SurveyA, SurveyB, delal=True, stochdrop=True):
  ''' abcpmc differs from astroABC in that sense that the model data is called before the data in the metric arguement,
      so it is metric(model,data) instead of metric(data,model)
      
  So this is just a wrapper    
  '''
  return ABC_dist_severalMetrices( SurveyB, SurveyA, delal=True, stochdrop=stochdrop)