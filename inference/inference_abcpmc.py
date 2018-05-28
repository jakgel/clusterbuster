#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 11:51:56 2017

@author: jakobg
"""
from __future__ import division, print_function

# start by importing astroabc and numpy
import numpy as np
import abcpmc
import corner
import time
import matplotlib.pyplot                   as plt
import clusterbuster.iout.misc             as iom
import surveysim.music2.runsurvey          as music2run
import inference.surveymetrics             as surmet
import seaborn as sns



''' The Test Part: Taken from http://nbviewer.jupyter.org/github/jakeret/abcpmc/blob/master/notebooks/dual_abc_pmc.ipynb'''

def create_new_sample(theta, samples_size=1000):
    mu,sigma = theta
    if sigma<=0:
        sigma=10
    return np.random.normal(mu, sigma, samples_size)


def dist_measure(x, y):
    return [np.abs(np.mean(x, axis=0) - np.mean(y, axis=0)),
            np.abs(np.var(x, axis=0) - np.var(y, axis=0))]


def main_test(Nthreads=7):
    
    
    import matplotlib
    import pandas as pd
    

    sns.set_style("white")
    np.random.seed(10)
    matplotlib.rc("font", size=30)

    
    samples_size = 1500
    N_varySame   = 2000
    mean = 2
    sigma = 1
    data = np.random.normal(mean, sigma, samples_size)
    print(data)
    
    ''' Plot 1 '''
    f,ax = plt.subplots()
    sns.distplot(data)
    
    
    ''' Plot 2 '''
    distances = [dist_measure(data, create_new_sample((mean, sigma), samples_size)) for _ in range(N_varySame)]
    dist_labels = ["mean", "var"]
    
    f, ax = plt.subplots()
    sns.distplot(np.asarray(distances)[:, 0], axlabel="distances", label=dist_labels[0])
    sns.distplot(np.asarray(distances)[:, 1], axlabel="distances", label=dist_labels[1])
    plt.title("Variablility of distance from simulations")
    plt.legend()
    
    
    ''' The abcpmc part starts '''
    prior = abcpmc.GaussianPrior(mu=[3.5, 2.0], sigma=np.eye(2) * 0.5)
    
    alpha = 85  
    T = 20
    eps_start = [1.0, 1.0]
    
    sampler = abcpmc.Sampler(N=1000, Y=data, postfn=create_new_sample, dist=dist_measure, threads=Nthreads)
    sampler.particle_proposal_cls = abcpmc.OLCMParticleProposal
    eps = abcpmc.ConstEps(T, eps_start)
    

    t0 = time.time()
    pools = launch(sampler,prior,alpha,eps)
    print("took", (time.time() - t0))
    
     


''' My own testpart '''
def testrand(notused, randomseed=False):
    if randomseed:
        import random
        seed = random.randrange(4294967295)
        np.random.seed(seed=seed)      
        print("Seed was:", seed)
#    print(np.random.poisson(4,5))

    return 0

def testmetric(inp1, inp2):
    ''' The testmetric does not care about the inputs 
        Returns: 3 random gaussian values
    '''

#    print(np.random.normal(0.5,0.1,3))
    return np.random.normal(1.0,0.15,3)




''' Functionalities '''

def plot_abctraces(pools, surveypath=''):
    
    ''' Input: a list of pools in the abc format 
        Generates trace plots of the thetas,eps and metrics '''



    ''' Plot Metric-Distances '''
    distances = np.array([pool.dists for pool in pools])
    print(distances.shape)
    f, ax = plt.subplots()
    for ii in  range(distances.shape[2]):
        ax.errorbar(np.arange(len(distances)), np.mean(distances, axis=1)[:, ii], np.std(distances, axis=1)[:, ii], label='M%i' % (ii+1))
#            sns.distplot(np.asarray(distances)[:, ii], axlabel="distances", label='M%i' % (ii))   
    plt.title("Development of Metric Distances")
    plt.xlabel('Iteration')
    plt.ylabel('Distance')
    plt.legend()
    plt.savefig('%s/Metrics.png' % (surveypath))

    ''' Plot Variables '''
    thetas = np.array([pool.thetas for pool in pools])
    print(thetas.shape)
    f, ax = plt.subplots()
    for ii in  range(thetas.shape[2]):
        ax.errorbar(np.arange(len(thetas)), np.mean(thetas, axis=1)[:, ii], np.std(thetas, axis=1)[:, ii], label='theta %i' % (ii+1))
#            sns.distplot(np.asarray(distances)[:, ii], axlabel="distances", label='M%i' % (ii))   
    plt.title("Development of Parameters")
    plt.xlabel('Iteration')
    plt.ylabel('Theta')
    plt.legend()
    plt.savefig('%s/Thetas.png' % (surveypath))
    
    ''' Plot Variables '''
    thetas = pool.thetas
    figure = corner.corner(thetas)
    figure.savefig('%s/CornerThetas_%02i.png' % (surveypath,len(pools)-1))
    
    '''
    corner.corner(distances)
    plots the various distances over each other and shows nicely that they are uncorrelated. 
    This is not super important, you could also use correlated distances with this abbroach. On the other hand it is interesting to see
    that both measures are independent, often this is a sign that they are good features!
    '''
    

    ''' Plot Epsilon'''
    fig,ax = plt.subplots()
    eps_values = np.array([pool.eps for pool in pools])
    for ii in  range(distances.shape[2]):
        ax.plot(eps_values[:, ii], label='M%i' % (ii))
    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"$\epsilon$", fontsize=15)
    ax.legend(loc="best")
    ax.set_title("Thresholds")
    fig.savefig('%s/Thresholds.png' % (surveypath))


        
    ''' Plot Parameters'''
    for ii,_ in enumerate(pools):
        jg = sns.jointplot(pools[ii].thetas[:, 0], 
                      pools[ii].thetas[:, 1], 
                      kind="kde", 
                     )
        jg.ax_joint.set_xlabel('var1')
        jg.ax_joint.set_ylabel('var2')
        jg.savefig('%s/FirstThetas_%i.png' % (surveypath, ii))

    return 0


def launch(sampler,prior,alpha,eps,ratio_min = 4e-2,surveypath=None, pool=None):
    ''' Launches '''
    ''' Could become impleneted in abcpmc '''
    
    pools = []
    for pool in sampler.sample(prior, eps, pool):
        eps_str = ", ".join(["{0:>.4f}".format(e) for e in pool.eps])
        print("T: {0}, eps: [{1}], ratio: {2:>.4f}".format(pool.t, eps_str, pool.ratio))

        for i, (mean, std) in enumerate(zip(*abcpmc.weighted_avg_and_std(pool.thetas, pool.ws, axis=0))):
            print(u"    theta[{0}]: {1:>.4f} \u00B1 {2:>.4f}".format(i, mean,std))

        eps.eps = np.percentile(pool.dists, alpha, axis=0) # reduce eps value
        pools.append(pool)
        
        iom.pickleObject(pools, surveypath, 'launch_pools', append = False)
        
        
        
        ''' Creates plots on the fly '''
        plot_abctraces(pools,surveypath)
 
        if pool.ratio < ratio_min:
            print('Ended abcpmc because ratio_min < %.3e' % (ratio_min))
            break
    sampler.close()
    return pools
    
    
    
    
'''=== Here starts the flesh and meet of the analysis done for my PhD ==='''

def main_MUSIC_2(dataNVSS, Nthreads=15, new_run=True):

    ''' 
    eps_start is actually important as the next iteration will only start if the number of computed trials within these boundaries will be Nw.
    So in one case I had to draw and compute twice as many particles than Nw.
    
    
    
    
    About the treads: 14-16 treads are fine, as more treats wont be fully used and just stuff the taskqueue
    
    Test & Development History: 
      05855 range(192):   16? iterations: ...         was done with screen. Terminated early ..., ...  
      10295 range(192):   > 1 iterations: ...         acceptance rate (0.80 instead of 0.875  & introduced maxtsksperchild & chunksize BUT should work nice!
      28914 range(240):   16? iterations: ... smaller acceptance rate (0.7 instead of 0.85 --> faster prunning and more noise), Flat prior. Ran into swap  before 80 computations
      21631 range(128):   12+ iterations ... a mix of models, a try towards a real model inference. Somehow worked (even with some memory overrun and no temrination like exspected)
      22326 range(40):  failed iteration ... amix of different models directly from the prior
      XXXX  range(20):      no iteration ... kappa    much  flatter than baseline
      27581 range(20):      no iteration ... kappa a little flatter than baseline
      9229  range (20)      no iteration ... the baseline model
    '''
    import matplotlib
    import pandas as pd
    surveypath = '/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS/'

    with open('/data/ClusterBuster-Output/AllSurveys/count.txt', 'w') as f:
        f.write('0')
    
    sns.set_style("white")
    #np.random.seed(10)
    matplotlib.rc("font", size=30)

    data      = dataNVSS
    alpha     = 75     # acceptance percentile i.e. best 85 percent
    T         = 12     # Maximum number of iterations
    Nw        = 135    # Number of walkers per iteration
    maxtasksperchild=1 #int(Nw/Nthreads/2)  # e.g. 252/14/2 = 9
    eps_start = [3, 1, 4, 2] 
    
    
    ''' The abcpmc part starts: Define thetas i.e. parameter values to be inferred '''
#    mu_lgeff, mu_lgB0, mu_kappa  = -5,0 , 0.0, 0.5  #-4.5,0.5 with full sky coverage
#    means    = [mu_lgeff, mu_lgB0, mu_kappa]
#    COV      = [[0.8, -0.5, 0], [-0.5, 0.8, 0], [0, 0, 0.3]]

#    prior = abcpmc.GaussianPrior(mu=means, sigma=COV)
    prior = abcpmc.TophatPrior([-6.5,-1.5,-1], [-4.5,2,1.5])
    

    '''PhD plots'''
#    eps_start =  [5.3, 5.3, 5.3 ]
    T                = 20
    maxtasksperchild = 3
    prior            = abcpmc.TophatPrior([-7], [-4])
    eps_start        = [5]            
#    mu_lgeff, mu_lgB0, mu_kappa  = -5.0, 0.0, 0.5
#    means    = [mu_lgeff, mu_lgB0, mu_kappa]
#    COV      = [[0.0, -0.0, 0], [-0.0, 0.0, 0], [0, 0, 0.0]]
#    prior = abcpmc.GaussianPrior(mu=means, sigma=COV)
    '''PhD plots' END'''
 
    
    eps = abcpmc.ConstEps(T, eps_start)
#    sampler = abcpmc.Sampler(N=Nw, Y=data, postfn=testrand, dist=testmetric, threads=Nthreads, maxtasksperchild=maxtasksperchild)
    sampler = abcpmc.Sampler(N=Nw, Y=data, postfn=music2run.main_ABC, dist=surmet.abcpmc_dist_severalMetrices, threads=Nthreads, maxtasksperchild=maxtasksperchild)
    sampler.particle_proposal_cls = abcpmc.OLCMParticleProposal
    
    ''' compare with AstroABC
    sampler = astroabc.ABC_class(Ndim,walkers,data,tlevels,niter,priors,**prop)
    sampler.sample(music2run.main_astroABC)    
    '''
    
    
    t0 = time.time()
    
    #	startfrom=iom.unpickleObject('/data/ClusterBuster-Output/MUSIC_NVSS02_Test01/launch_pools')
    pool     = None #startfrom[-1]
    pools = launch(sampler,prior,alpha,eps,surveypath=surveypath,pool=pool)
    print("took", (time.time() - t0))


    return 0

if __name__ == "__main__":
    folder    = '/data/ClusterBuster-Output/'
    dataNVSS   = iom.unpickleObject('%s%s/pickled/Survey' % (folder,'NVSS'))
#    dataMUSIC2 = iom.unpickleObject('%s%s/pickled/Survey' % (folder,'MUSIC2COOL_NVSS_SSD_00002'))
#    surmet.abcpmc_dist_severalMetrices( dataMUSIC2, dataNVSS, delal=False)
    main_MUSIC_2(dataNVSS, new_run=True, Nthreads=20)
