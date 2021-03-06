#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 11:51:56 2017

@author: jakobg
"""
from __future__ import division, print_function

# start by importing astroabc and numpy
import json
import argparse
import configparser

import numpy as np
import abcpmc
import corner
import matplotlib
import matplotlib.pyplot          as plt
import clusterbuster.iout.misc    as iom
import surveysim.music2.runsurvey as music2run
import inference.surveymetrics    as surmet
import seaborn as sns

from functools import partial

"""=== Here starts the flesh and meet of the analysis done for my PhD ==="""
def main_abcpmc_MUSIC2(conf, test=False):

    """ 
    config should contain
    [][]: a list etc.
    
    
    eps_start is actually important as the next iteration will only start if the number of computed trials within these boundaries will be Nw.
    So in one case I had to draw and compute twice as many particles than Nw.
    
    About the treads: 14-16 treads are fine, as more treats wont be fully used and just sit in the taskqueue
    """
    
    
    # Loads the real data to compare with (and if neccessary also test data)
    data = iom.unpickleObject(conf['paths']['surveyreal'])
    if test:
        dataMUSIC2 = iom.unpickleObject(conf['paths']['surveysim'])
        print(type(dataMUSIC2.Rmodel), conf['paths']['surveysim'])
        surmet.abcpmc_dist_severalMetrices(dataMUSIC2, data, metrics=json.loads(conf['metrics']['used']),
                                           delal=False, stochdrop=conf['flavor']['stochdrop'],
                                           phoenixdrop = conf['flavor']['phoenixdrop'], outpath='/data/')
        return 0

    """ The abcpmc part starts: Define thetas i.e. parameter values to be inferred  and priors"""

    if conf['prior']['type'] == 'tophat':
        bounds = json.loads(conf['prior']['bounds'])
        prior = abcpmc.TophatPrior(bounds[0], bounds[1])
    elif conf['prior']['type'] == 'gaussian':
        means = json.loads(conf['prior']['means'])
        COV = json.loads(conf['prior']['covariance'])
        prior = abcpmc.GaussianPrior(mu=means, sigma=COV)
    else:
        print('inference_abcpmc::main_abcpmc_MUSIC2: prior %s is unknown!' % (conf['prior']['type']))
        return 0

    eps = abcpmc.ConstEps(conf.getint('pmc', 'T'), json.loads(conf['metrics']['eps_startlimits']))

    if test:
        sampler = abcpmc.Sampler(N=conf.getint('pmc', 'Nw'), Y=data, postfn=testrand, dist=testmetric,
                                 threads=conf.getint('mp', 'Nthreads'), maxtasksperchild=conf.getint('mp', 'maxtasksperchild'))
    else:
        sampler = abcpmc.Sampler(N=conf.getint('pmc', 'Nw'), Y=data, postfn=partial(music2run.main_ABC, parfile=conf['simulation']['parfile']),
                                 dist=partial(surmet.abcpmc_dist_severalMetrices, metrics=json.loads(conf['metrics']['used']),
                                              outpath=conf['paths']['abcpmc'], stochdrop=conf['flavor']['stochdrop'],
                                              phoenixdrop = conf['flavor']['phoenixdrop']),
                                 threads=conf.getint('mp', 'Nthreads'), maxtasksperchild=conf.getint('mp', 'maxtasksperchild'))

        # Prepares the file for counting
        with open(conf['paths']['abcpmc'] + 'count.txt', 'w+') as f:
            f.write('0')

    sampler.particle_proposal_cls = abcpmc.OLCMParticleProposal
    
    """ compare with AstroABC
    sampler = astroabc.ABC_class(Ndim,walkers,data,tlevels,niter,priors,**prop)
    sampler.sample(music2run.main_astroABC)    
    """
    #	startfrom=iom.unpickleObject('/data/ClusterBuster-Output/MUSIC_NVSS02_Test01/launch_pools')
    pool = None #startfrom[-1]
    launch(sampler, prior, conf.getfloat('pmc','alpha'), eps, surveypath=conf['paths']['abcpmc'], pool=pool)


def launch(sampler, prior, alpha, eps, ratio_min=1e-2, surveypath=None, pool=None, plotting=False):
    """ Launches pools
     Could become implemented in abc-pmc itself"""
    
    pools = []
    for pool in sampler.sample(prior, eps, pool):
        eps_str = ", ".join(["{0:>.4f}".format(e) for e in pool.eps])
        print("T: {0}, eps: [{1}], ratio: {2:>.4f}".format(pool.t, eps_str, pool.ratio))

        for i, (mean, std) in enumerate(zip(*abcpmc.weighted_avg_and_std(pool.thetas, pool.ws, axis=0))):
            print("    theta[{0}]: {1:>.4f} +- {2:>.4f}".format(i, mean, std))

        eps.eps = np.percentile(pool.dists, alpha, axis=0) # reduce eps value
        pools.append(pool)
        
        iom.pickleObject(pools, surveypath, 'launch_pools', append=False)

        """ Creates plots on the fly """
        if plotting:
            plot_abctraces(pools, surveypath)
 
        if pool.ratio < ratio_min:
            print('Ended abc-pmc because ratio_min < %.3e' % (ratio_min))
            break
    sampler.close()
    return pools
    



""" My own testpart """
def testrand(notused, randomseed=False):
    if randomseed:
        import random
        seed = random.randrange(4294967295)
        np.random.seed(seed=seed)      
        print("Seed was:", seed)
#    print(np.random.poisson(4,5))

    return 0

def testmetric(inp1, inp2):
    """ The testmetric does not care about the inputs 
        Returns: 3 random gaussian values
    """

#    print(np.random.normal(0.5,0.1,3))
    return np.random.normal(1.0,0.15,3)




""" Functionalities """

def plot_abctraces(pools, surveypath=''):
    
    """ Input: a list of pools in the abc format 
        Generates trace plots of the thetas,eps and metrics """

    sns.set_style("white")
    matplotlib.rc("font", size=30)

    """ Plot Metric-Distances """
    distances = np.array([pool.dists for pool in pools])
    print(distances.shape)
    f, ax = plt.subplots()
    for ii in range(distances.shape[2]):
        ax.errorbar(np.arange(len(distances)), np.mean(distances, axis=1)[:, ii], np.std(distances, axis=1)[:, ii], label='$\\Delta_%i$' % (ii+1))
#            sns.distplot(np.asarray(distances)[:, ii], axlabel="distances", label='M%i' % (ii))   
    #plt.title("Development of Metric Distances")
    plt.xlabel('Iteration')
    plt.ylabel('Distance $\Delta$ in metric')
    plt.legend()
    plt.savefig('%s/Metrics.png' % (surveypath))

    """ Plot Variables """
    thetas = np.array([pool.thetas for pool in pools])
    print(thetas.shape)
    f, ax = plt.subplots()
    for ii in range(thetas.shape[2]):
        ax.errorbar(np.arange(len(thetas)), np.mean(thetas, axis=1)[:, ii], np.std(thetas, axis=1)[:, ii], label='$\\theta_%i$' % (ii+1))
    #plt.title("Development of Parameters")
    plt.xlabel('Iteration')
    plt.ylabel('$\\theta_i$')
    plt.legend()
    plt.savefig('%s/Thetas.png' % (surveypath))
    
    
    """ Plot Variables """
    #TODO: Fix bug ... you need to call pools or pool?
    for ii, pool, in enumerate(pools):
        thetas = pool.thetas
        figure = corner.corner(thetas)
        plt.savefig('%s/CornerThetas_%02i.png' % (surveypath,ii))
    
    """
    corner.corner(distances)
    plots the various distances over each other to show if they are uncorrelated. 
    This is not super important, you could also use correlated distances with this approach. On the other hand it is interesting to see
    if both metrices are independent, often this is a sign that they are good features!
    """
    

    """ Plot Epsilon"""
    fig, ax = plt.subplots()
    eps_values = np.array([pool.eps for pool in pools])
    for ii in range(distances.shape[2]):
        ax.plot(eps_values[:, ii], label='$\epsilon_%i$' % (ii))
    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"$\epsilon$", fontsize=15)
    ax.legend(loc="best")
    #ax.set_title("Thresholds $\epsilon$")
    plt.savefig('%s/Thresholds.png' % (surveypath))


    """ Violin Plots """
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
    
    # generate some random test data
    all_data = [pool.thetas[:,0] for pool in pools]
    len_pools = [pool.thetas.shape[0]/pool.ratio for pool in pools]
    print('len_pools:', len_pools, sum(len_pools))
    mod_data = [np.concatenate((pool.thetas[:,0]+0.2, pool.thetas[:,0]), axis=0) for pool in pools]
    # plot violin plot
    #background = axes.violinplot(mod_data,
    #                   showmeans=False,
    #                   showmedians=False, showextrema=False)
    axes.violinplot(all_data,
                       showmeans=False,
                       showmedians=True)
    
    #for pc in background['bodies']:
    #    pc.set_facecolor('grey')
#        pc.set_edgecolor('black')
    #    pc.set_alpha(0.4)
        
    #for pc in foreground['bodies']:
    #    pc.set_facecolor('cornflowerblue')
#        pc.set_edgecolor('black')
    #    pc.set_alpha(1)

#    axes.set_title('violin plot')
    
    # adding horizontal grid lines
    axes.yaxis.grid(True)
    axes.set_xticks([y+1 for y in range(len(all_data))])
    axes.set_xlabel('Iteration')
    axes.set_ylabel('$\\log_{10}(\\xi_e)$')
    axes.set_ylabel('$\\log_{10}(\\xi_e)$')

    # add x-tick labels
    plt.setp(axes, xticks=[y+1 for y in range(len(all_data))])
    plt.savefig('%s/Violin.png' % (surveypath))
    plt.savefig('%s/Violin.pdf' % (surveypath))
    plt.clf()

    """ Plot Parameters
    pools[ii].thetas[:, 0] is a numpy array 
    """

    for ii, nouse in enumerate(pools):
        if thetas.shape[1] > 1:
            jg = sns.jointplot(pools[ii].thetas[:, 0], 
                          pools[ii].thetas[:, 1], 
                          #kind="kde",  # BUG: creates an error
                         )
            jg.ax_joint.set_xlabel('var1')
            jg.ax_joint.set_ylabel('var2')
            plt.savefig('%s/FirstThetas_%i.png' % (surveypath, ii))

    return 0


    


if __name__ == "__main__":
    import clusterbuster.iout.misc as iom
    #surveypath ='/data/ClusterBuster-Output/abcpmc-MUSIC2NVSS_Run_06/'
    #pools = iom.unpickleObject(surveypath+'launch_pools')
    #print('Acceptance ratios:', [pool.ratio for pool in pools])

    #plot_abctraces(pools, surveypath=surveypath)
    #exit()

    """ Full routines for parsing a combination of argparse, configparser and json"""
    parser = argparse.ArgumentParser(description='Evaluates the best solutions of survey simulations with the abc approach.')
    parser.add_argument('-parini', dest='parini', action='store', default='params.ini', type=str,
                        help='set filepath for parameter initialisation file.')
    args = parser.parse_args()

#    fd = open(args.parini, 'r')
    config = configparser.ConfigParser()
    config.read(args.parini)
    args = parser.parse_args()
    
    main_abcpmc_MUSIC2(config, test=False)  #dataNVSS
