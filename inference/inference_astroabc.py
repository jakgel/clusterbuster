#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 11:51:56 2017

@author: jakobg
"""

# start by importing astroabc and numpy
import numpy as np
import astroabc
import time
import matplotlib.pyplot                   as plt
import pyutil_my.IOutil                    as iout
import pyutil_my.SurveyMetrics             as SurveyMetrics
import Analysis_MUSIC2.RunSurvey           as MUSIC2

#define a method for simulating the data given input parameters
def model_sim(param, pool=None):
    #Ideally do something with the pool here
    cov =np.array([0.01,0.005,0.005,0.1])
    return astroabc.Model("normal",10000).make_mock((param,cov))

def model_sim_mod(param, pool=None):
    param =  [param[0][0],param[0][1]]
#    print param
    #Ideally do something with the pool here
    cov =np.array([0.01,0.005,0.005,0.1])
    return astroabc.Model("normal",10000).make_mock((param,cov))

def dist_metric(d,x):
    return np.sum(np.abs(np.mean(x,axis=0) - np.mean(d,axis=0)))


#=====
def main_test():
    import plotly.plotly as py
    import plotly.graph_objs as go
    
    #make the fake data with diagonal covariance
    means= np.array([0.037579, 0.573537])
    cov =np.array([0.01,0.005,0.005,0.1])
    COV = [[0.05, 0.025], [0.025, 0.5]]
#    data = astroabc.Model("normal",100).make_mock((means,cov))
    data = astroabc.Model("normal",100).make_mock(([1,1],cov))
    case = 3




    """ Implemented """
    if case == 1:
        priors = [('uniform', [-0.1,0.2]), ('uniform', [0.1, 0.9])]
        Nparas = 2

    """Issue 5: var A; not working """
    if case == 2:
        
        multivariate = ('multigauss', [means[0],means[1],COV])
        priors = [multivariate]
        Nparas = 1
        
    """Issue 5: var B; working wekk """
    if case == 3:
        from scipy.stats import multivariate_normal
        
        prior =np.random.multivariate_normal(means, COV, 200) #cov.reshape((2,2))
        np.savetxt('output_1.txt', prior, fmt='%+0.4e')
        priors = [('nonstandard', ["output_",1,0]),('nonstandard', ["output_",1,1])]  #('uniform', [value_1,value_2]), 
        Nparas = 2


        """For the non-standard option the format expected is from a single pdf for one parameter written as text column output. 
           Let's say you have the output from an MCMC run, each column in the output file would then correspond to the marginalized 1D pdf for each of 
           the parameters. E.g. for parameter A whose pdf is in column 0 in "output_1.txt" use:
                          
                          NFILES = 1;
                          FILENAME = "output_"    
                          
                          and prior setting priors = [('nonstandard', [1])] 
        """

                                    

    prop   = {'dfunc':dist_metric, 'outfile':"gaussian_exampleB.txt", 'verbose':1, 
              'adapt_t': True, 'variance_method':4, 'k_near':10, 'mp': True,  'num_proc':20 }

    walkers = 160
    steps   = 5
    sampler = astroabc.ABC_class(Nparas,walkers,data,[0.5,0.002],steps,priors,**prop)

    if case == 2:
        sampler.sample(model_sim_mod)
    else:
        sampler.sample(model_sim)

    A = np.loadtxt('gaussian_exampleB.txt',skiprows=1)
    plot_chain2D(prior, A, walkers, steps, extent = (-0.5,0.5,-1.0,2.0)) 
    
    return 0

""" From https://plot.ly/pandas/contour-plots/ """
def kde_scipy( m1, m2, extent=None ):
    import scipy.stats as st
    # http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.gaussian_kde.html
    if extent is None:
        xmin = m1.min()
        xmax = m1.max()
        ymin = m2.min()
        ymax = m2.max()
    else:
        xmin, xmax, ymin, ymax = extent

    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = st.gaussian_kde(values)
    #print X
    Z = np.reshape(kernel(positions).T, X.shape).T

    return Z

def plot_chain2D (prior, A, walkers, steps, extent = (-0.5,0.5,-1.0,2.0), outfile=None):
    import corner
 # This part somehow is buggy ... if I shift the y-axis the x-axis is shifted but the data is ok
 
    figure = corner.corner(prior, quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 12})
 
 
    plt.figure()
    plt.title('Step %i in the MCMC-ABC chain (prior)' % (0))
        
    X = prior[:,0]
    Y = prior[:,1]
    Z = kde_scipy( X, Y, extent = extent)
    
    im = plt.imshow(Z, interpolation='bilinear', origin='lower',extent=extent, aspect='auto')
    plt.colorbar()
    CS = plt.contour(Z, origin='lower',
                         linewidths=2,extent=extent)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel('$\lg f_e$')
    plt.ylabel('$\lg B_0$')
    
    for ii,count in enumerate(np.arange(steps)):
        X = A[(count)*walkers:(count+1)*walkers,0]
        Y = A[(count)*walkers:(count+1)*walkers,1]
        Z = kde_scipy( X, Y, extent = extent )

        plt.figure()
        plt.title('Step %i in the MCMC-ABC chain' % (ii+1))
        
        im = plt.imshow(Z, interpolation='bilinear', origin='lower',extent=extent, aspect='auto')
        plt.colorbar()
        CS = plt.contour(Z, origin='lower',
                         linewidths=2,extent=extent)
        plt.clabel(CS, inline=1, fontsize=10)
        plt.xlabel('$\lg f_e$')
        plt.ylabel('$\lg B_0$')
    return 0


def plot_chain2D_new(prior, A, walkers, steps, extent = (-0.5,0.5,-1.0,2.0), outfile=None):
 

    X = prior[:,0]
    Y = prior[:,1]
    Z = kde_scipy( X, Y, extent = extent)
    

    import pylab as pyl
    
    NColumns = 3
    number_of_subplots= steps+1 + ((steps+1)%NColumns)

    # alternatively https://stackoverflow.com/questions/20057260/how-to-remove-gaps-between-subplots-in-matplotlib
    
    ax0 = pyl.subplot(number_of_subplots,NColumns,1)
    im = plt.imshow(Z, interpolation='bilinear', origin='lower',extent=extent, aspect='auto') 
    CS = ax0.contour(Z, origin='lower',linewidths=2,extent=extent)
    ax0.xaxis.set_visible(number_of_subplots-(1)<NColumns)
    ax0.text(0.97, 0.72, 'Prior',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax0.transAxes,
        color='w', fontsize=6)
    
    for ii,count in enumerate(np.arange(steps)): #xrange(number_of_subplots)):
        X = A[(count)*walkers:(count+1)*walkers,0]
        Y = A[(count)*walkers:(count+1)*walkers,1]
        Z = kde_scipy( X, Y, extent = extent )

        ax1 = pyl.subplot(number_of_subplots,NColumns,count+2, sharey=ax0, sharex=ax0)
        im = ax1.imshow(Z, interpolation='bilinear', origin='lower',extent=extent,aspect='auto') #, aspect='auto'
        ax1.xaxis.set_visible(number_of_subplots-(count+2)<=NColumns)
        ax1.yaxis.set_visible((count+1)%NColumns==0)
        CS = ax1.contour(Z, origin='lower',linewidths=2,extent=extent)
        ax1.clabel(CS, inline=1, fontsize=4)
        ax1.text(0.97, 0.72, '%i' % (count+1),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax1.transAxes,
            color='w', fontsize=6)
    cax = plt.axes([0.92, 0.58, 0.035, 0.3])  
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.colorbar(im, cax=cax)
    plt.savefig('test.pdf')
    plt.show()

    return 0





      

def main_MUSIC_2(newrun=True):
    """ Runs astroABC on MUSIC-2 vs NVSS"""
    folder     = '/data/ClusterBuster-Output/'
    
    


    start = time.time()

    if newrun:
        walkers  = 100 #100       # number of particles/walkers
        niter    = 5 #4           # number of iterations
        tlevels  = [3.0,0.9]      # maximum,minimum tolerance
        Ndim     = 3
        mu_lgeff, mu_lgB0  = -5.0 , -0.5  #-6.0,-0.5 with less sky coverage
        means    = [mu_lgeff, mu_lgB0]
        COV      = [[0.8, -0.5],[-0.5, 0.8]]



#        """ SuperMegaTest """
#        mu_lgeff, mu_lgB0  = -5.0 , -0.5  #-6.0,-0.5 with less sky coverage
#        COV      = [[0.001, -0.0], [-0.0, 0.001]] 
#        walkers  = 16
#        niter    = 1

        """ Caseof a non-standart-prior """
    
        prior =np.random.multivariate_normal(means, COV, 20000) #cov.reshape((2,2))
        np.savetxt('output_1.txt', prior, fmt='%+0.4e')
        priors = [('nonstandard', ["output_",1,0]),('nonstandard', ["output_",1,1]), ("uniform",[0,1.0])]  #('uniform', [value_1,value_2]), 
    
    
        """  This is the pre-existing electron part Gelszinnis_Model """
#        N_pre    = 4 
#        Ndim    += N_pre
#        priors  += [("normal",[-4,2]),("uniform",[0,1.0]),("normal",[-4,2]),("uniform",[0.3,2.0])] #lgp0, p_sigma, lgsigmoid_0, sigmoid_width)
        outfile = folder+"abc_pmc_output_run_05_"+str(Ndim)+"param.txt"
        
                                                     
        """  This is the pre-existing electron part Hoeft_Model 
        
        I had to made a decicion about some values, so I tried to fix ...
        B0,kappa,eff
        t0,t1    = 0.1,5 --> varying   ---> get these values as logarithmic values? better not
        Also fix t0 as 10% of t1, because if not, you could get an odd correlation ...
        
        ratio    = 0.05 --> varying
        nu1,nu2  = 1e-6,1e-2
        
        
        Metric: Add mach number and maybe the disance of the mean of the othe rcorrelations?
           ---> Somehow like the entropy of the correlation matrix ? ...
           ---> Prewhiten The data (zero, mean, standard etc. for NVSS)
           ---> Transform the simulated data in the saem frame and measure the distance**2 along PCA1,PCA2,PCA3 as metrices etc, ...
           --> This wont work, because data with zero scatter would be good. 
           --> Better take the 
        """

    #
        prop={'tol_type':'exp',"verbose":1,'adapt_t':True,
			 'threshold':75, 'variance_method':4, 'k_near':10,    #'threshold':75,'pert_kernel':2,'variance_method':0,
              'dist_type': 'user','dfunc':SurveyMetrics.ABC_dist_music2, 'restart':folder+"restart_test.txt", 
              'outfile':outfile,'mpi':False,
              'mp':True,'num_proc':16, 'from_restart':False}
        
        """Loads NVSS survey"""
        surveyname = 'NVSS'
        data = iout.unpickleObject('%s%s/pickled/SurveySample' % (folder,surveyname))
             
        

        print priors
    
    
    #    multivariate = ('multigauss', [mu_lgeff,mu_lgB0,COV])
    #    prior = [multivariate, ("uniform",[0,1.0])]
                  
        sampler = astroabc.ABC_class(Ndim,walkers,data,tlevels,niter,priors,**prop)
        sampler.sample(MUSIC2.main_astroABC)      
    


    #oldrun
    else:
        
        run = 4
        means    = [-5.0, -0.5]
        COV      = [[0.8, -0.5], [-0.5, 0.8]]
        
        
        """ run 7: First try with (fake) preexisting pasma
            + From here on MUSIC2_NVSS02_SSD.parset (with compression enhanced B-field) """
        if run == 7:
            tlevels  = [3.0,0.6]
            means    = [-5.0,-0.5]
            walkers  = 100
            niter    = 6
            outfile = folder + 'abc_pmc_output_run_06_7param.txt'
            
            """
            standard metric (pls copy here)
            """
            
        """ run 6: Test only one (best) model without any covariance """
        if run == 5:
            tlevels  = [3.0,0.9]
            means    = [-5.0,-0.5]
            COV      = [[0.001, -0.0], [-0.0, 0.001]] 
            walkers  = 20
            niter    = 0
            outfile = folder + 'abc_pmc_output_run_05_3param.txt'
            
            """
            standard metric (pls copy here)
            """     
            
        """ run 5: From here on MUSIC2_NVSS02_SSD.parset(4 times larger, more particles with lower mach numbers)
                   From here in with compression enhanced B-field"""
        if run == 5:
            tlevels  = [3.0,0.9]
            means    = [-5.0,-0.5]
            walkers  = 100
            niter    = 5
            outfile = folder + 'abc_pmc_output_run_05_3param.txt'
            
            """
            standard metric (pls copy here)
            """    
            
        
        """ run 4 """
        if run == 4:
            tlevels  = [3.0,0.9]
            means    = [-5.0,0.0]
            walkers  = 100
            niter    = 9
            outfile = folder + 'abc_pmc_output_run_04_3param.txt'
            
            """
            standard metric (pls copy here)
            """
        
        """ run 3 """
        if run == 3:
            tlevels  = [3.0,0.6]
            means    = [-5.0,0.0]
            walkers  = 300
            niter    = 11
            outfile = folder + 'abc_pmc_output_run_03_3param.txt'
                
            """
            standard metric (pls copy here)
            did terminate at step 5/6 because no of the samples fulfiled the strict rejection criterium
            """
        
        """ run 2: From here in survey completeness 0.85 ?"""
        if run == 2:
            tlevels = [3.0,0.9]
            means   = [-5.0,0.0]
            walkers = 180
            niter   = 8
            outfile = folder +  'abc_pmc_output_run_02_3param.txt'
            """standard metric (pls copy here)""" 


        """ run 1: survey completeness 0.25"""
        if run == 1:
            outfile = folder + 'abc_pmc_output_3param_0.25simple.txt'  
            tlevels  = [3.0,1.1]      # maximum,minimum tolerance      
    #        prior    = np.random.multivariate_normal(means, COV, 20000) #cov.reshape((2,2))
        A = np.loadtxt(outfile,skiprows=1)
        print len(A[:,0:2])      
    #        plot_chain2D(prior, A[:,0:2], 70, 6, extent=(-5.5,-3,-1.2,1.7), outfile=outfile+'.pdf')
    #        plot_chain2D_new(prior, A[:,0:2], 70, 6, extent=(-5.5,-3,-1.2,1.7), outfile=outfile+'.pdf')
        
        import corner    
        figure = corner.corner(A[:,0:4], labels=[r"$\lg f_e$", r"$lg B_0$", r"$\kappa$", r"Dist$_\mathrm{metric}$"],
                       quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
        plt.savefig(folder+'corner_allruns')
        plt.figure()

        
    end = time.time()
    print  'Elapsed time for main_MUSIC_2(): %.2e seconds.' %  (end - start)
    
    
    return 0



if __name__ == "__main__":
#    folder    = '/data/ClusterBuster-Output/'
#    dataNVSS  = iout.unpickleObject('%s%s/pickled/SurveySample' % (folder,'NVSS'))
#    dataMUSIC = iout.unpickleObject('%s%s/pickled/SurveySample' % (folder,'MUSIC2_NVSS01_14168'))  
#    ABC_dist_music2( dataNVSS, dataMUSIC)
    main_MUSIC_2(newrun=True)
#    main_test()
#