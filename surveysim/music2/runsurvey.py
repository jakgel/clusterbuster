#!/usr/bin/env python

"""
Until the 18 Apr 2017 there is no version of abcpmc that allows for np.random seeds within multiprocessing processes! 
The alternative astroABC that can handle this. One simple workaround is to set the np.seed(...) within the function that has to be called.
"""

from __future__ import division,print_function
import argparse
import time
import datetime
import os
import copy
import traceback

""" In the local folder ... """
import surveysim.music2.mockobs    as mockobs
import surveysim.music2.loadsnap   as loadsnap
import surveysim.music2.radiomodel as radiomodel   
""" ==="""

import multiprocessing as mupro
import numpy           as np
import pandas          as pd

import clusterbuster.surveyclasses  as cbclass
import clusterbuster.iout.misc      as iom
import clusterbuster.surveyut       as suut
import clusterbuster.constants      as myu


import time
import timeit
import random 
from random import uniform
from astropy.cosmology import FlatLambdaCDM



#==== 
par = argparse.ArgumentParser()
par.add_argument('--par',                 default='MUSIC2COOL_NVSS',  help='parset file to run' )
""" "--par 'MUSIC_NVSS01.parset" """

args, unknown = par.parse_known_args()

#=========================================================================================================================================================  
POISON_PILL = "STOP"


"""Originally main() was designed as independent functions because I wanted to prevent the buffer to overflow.
   Now I use these two functions as the one that controls the ABC (main) and the one that does the computation (shedule task)"""

def main(parfile, workdir=None, ABC=None, verbose=False, survey=None, index=None, Clfile='clusterCSV/MUSIC2-AGN',
         processTasks=True):

    """
    input: parameter file:
           ABC  : None or a list of parameters (what a pity that a dictionary is not possible with the current abcpmc or ABCpmc module)

    output: some strings to track down the output directories

    once a survey is not None all the default survey values will be used;
    i.e. the parfile will be considered much
    For future: - delete survey findings (or add for different detection parameter)
                or make the Rmodel changeable

    If pase['default'] is False then a run will be made with one shell and one single cluster realisation per snapshot

    For abcpmc MUSIC-2 vanilla use MUSIC2-AGN for clusterCSV and MUSIC2_NVSS02_SSD.parset for parset
               MUSIC-2cooling  use MUSIC2-AGN for clusterCSV and MUSIC2COOL_NVSS.parset         for parset
    """
    RModelID = os.getpid()  # get process id

    seed = random.randrange(4294967295)
    np.random.seed(seed=seed)

    """also possible : np.random.seed(seed=None) 
    ... this here is only to save the seed, which is not needed, because all the random stuff can be reconstructed from the cluster statistics

    processes from random are well seeded even without this, somehow numpy needs this kind of seed
    but only in abcpmc and not ABC

    import random
    random.randrange(sys.maxsize)
    """

    # ===  Read parset; then extract fundamental parameters ... Parset OBLIGATORY has be saved as  'parsets/*.parset'
    if workdir is None: workdir = os.path.dirname(__file__) + '/'
    pase, Z = suut.interpret_parset(parfile, repository=workdir + '/parsets/')
    # TODO: workdir should be the same as pase['miscdir'] yet the issue is that it is unknown before miscdir is known.

    if survey is None:
        surveyN = parfile.replace('.parset', '')
        savefolder = pase['outf'] + surveyN
        logfolder = pase['outf'] + surveyN
    else:
        surveyN = survey.name
        savefolder = survey.outfolder
        logfolder = survey.logfolder

    # Constants
    # ==== Cosmological Parameter - MUSIC-2
    Omega_M = 0.27  # Matter density parameter

    # === Create folder if needed
    iom.check_mkdir(savefolder)
    iom.check_mkdir(logfolder)

    smt = iom.SmartTiming(rate=pase['smarttime_sub'], logf=savefolder + '/smt');
    smt(task='Prepare_Clusters')

    # === Create detection information
    dinfo = cbclass.DetInfo(beam=[float(pase['S_beam']), float(pase['S_beam']), 0], spixel=float(pase['S_pixel']),
                            rms=float(pase['RMSnoise']) * 1e-6,
                            limit=float(pase['RMSnoise']) * float(pase['Detthresh']) * 1e-6,
                            nucen=float(pase['nu_obs']), center=(0, 0), survey='UVcoverage')

    if survey is None:
        """Create all galaxy clusters:
            All these steps are to decide which clusters to use. Better placement: in SurveyUtils"""

        """ Read the cluster lists from MUSIC-2 for all snapshots """
        all_clusters = pd.read_csv('%s%s_allclusters.csv' % (pase['miscdata'], Clfile))
        if verbose:
            all_clusters.info()
        zsnap_list = pd.Series(all_clusters['redshift'].unique())
        snapidlist = pd.Series(all_clusters['snapID'].unique())
        clusterIDs = list(all_clusters['clID'].unique())
        NclusterIDs = [len(all_clusters[all_clusters['redshift'] == z]) for z in zsnap_list]
        misslist = np.loadtxt(pase['miscdata'] + pase['missFile'])
        """ e.g. cluster 10# was not (re)simulated, in neither of the MUSIC-2 simulations (also 7 has some issues?)"""

        GClList = []

        if pase['snaplistz'] != 'None':

            snaplistz = [float(z) for z in iom.str2list(pase['snaplistz'])]
            snapidlist = [float(z) for z in iom.str2list(pase['snapidlist'].replace(' ', ''))]

            zsnap_list = (snaplistz[::-1])[0:11]  # 0:17 List of available sn [0:11]
            snapidlist = (snapidlist[::-1])[0:11]  # 0:17 [0:11]

            NclusterIDs = [len(all_clusters['clID'].unique().tolist())] * len(snaplistz)
            if verbose: print('NclusterIDs', NclusterIDs[0])

        use_list = [True] * len(zsnap_list)  # Also put some False, you don't really want to use the z=4.0 snapshots!
        Vsimu = (1.0 / (myu.H0 / 100.)) ** 3  # Gpc**3 comoving volume

        """ Iterate trough each shell of your z-onion and attribute clusters to them
            with z-range and percentage of covered sky, we have z=0.1
        """
        if not suut.TestPar(pase['default']):
            N_shells = 1
        else:
            N_shells = float(pase['N_shells'])

        shells_z, delta_z = np.linspace(Z[0], Z[1], num=N_shells + 1, retstep=True)
        cosmo = FlatLambdaCDM(H0=myu.H0, Om0=Omega_M)
        DCMRs = cosmo.comoving_volume(shells_z).value / 1e9
        count = 0

        for (zlow, zhigh, VCMlow, VCMhigh) in zip(shells_z[0:-1], shells_z[1:], DCMRs[0:-1], DCMRs[1:]):

            """ Iterate through each shell of the observed volume 
                and assign clusters
            """

            boundaries_z = (zlow, zhigh)
            VCM = (VCMlow, VCMhigh)
            z_central = np.mean(boundaries_z)
            if not suut.TestPar(pase['default']):
                z_central = 0.051
            choosen = suut.assign_snaps(zsnap_list, boundaries_z, VCM[1] - VCM[0], NclusterIDs,
                                        sigma_z=float(pase['sigma_z']),
                                        skycoverage=float(pase['surv_compl']), Vsimu=Vsimu, use_list=use_list,
                                        fake=(not suut.TestPar(pase['default'])), logmode=None)

            for (snap, kk) in choosen:
                l = all_clusters[(all_clusters["clID"] == clusterIDs[kk])
                                 & (all_clusters["snapID"] == snapidlist[snap])]
                """ It would be good if you could directly access the element like in an ordered list, 
                    as this would dramatically speed up the process
                """

                """ Skips missing snapshots --> they will also miss in the .csv"""
                if len(l) == 0:
                    if verbose:
                        print('__ Missing snapshot:', clusterIDs[kk], snapidlist[snap])
                    continue

                ids = int(l["clID"])
                M200 = float(l["M200"])

                # Filter for the cluster masses. Please mind that this filtering step is also redshift dependent
                if suut.TestPar(pase['empicut']) and np.log10(M200) < (13.6 + 2 * z_central):
                    continue

                count += 1

                # Decide on the projection of the cluster
                # it would be great if a random initializer between 0 and 1 could have been saved,
                if suut.TestPar(pase['rotation']):
                    theta = np.arccos(uniform(0, 2) - 1)
                    phi = uniform(0, 2 * np.pi)
                    psi = uniform(0, 2 * np.pi)
                else:
                    theta = 0
                    phi = 0
                    psi = 0

                # Create mockObs and the galaxyCluster_simulation
                mockObs = cbclass.MockObs(count, theta=theta, phi=phi, psi=psi, snap=snapidlist[snap],
                                          z_snap=zsnap_list[snap], clid=ids,
                                          snapfolder=pase['xrayfolder'], xrayfolder=pase['xrayfolder'],
                                          headerc=pase['headerC'])
                GClList.append(
                    cbclass.Galaxycluster_simulation("MUSIC2%05i-%06i-%06i" % (ids, snapidlist[snap], count), count,
                                                     z=z_central, M200=M200, dinfo=dinfo,
                                                     mockobs=mockObs))  # , **addargs

        # Also create a list of the chosen clusters for later loockup
        GClList = sorted(GClList, key=lambda x: (x.mockobs.clid, -x.mockobs.snap))

        if verbose: print('Length of GClList:', len(GClList))
        """ New approach: Create a list of modes for the radio emission (Rmodels)
        """
        surmodel = None
        if ABC is None:
            """CAVEAT: The issue with the currently used ABC routines is that you have to give them arrays. Which is why 
               the corresponding model associated has to be defined at this layer.
               Giving the procedure a function which would create this array would allow all of this to be defined 
               in the top layer of ABC
            """

            RModel = cbclass.RModel(RModelID, effList=[float(pase['eff'])], B0=float(pase['B0']),
                                    kappa=float(pase['kappa']), compress=float(pase['compress']))

            if suut.TestPar(pase['redSnap']):
                RModel.effList = RModel.effList[0]
        elif len(ABC) == 1:
            """ Vary only efficiency """
            (lgeff) = ABC
            RModel = cbclass.RModel(RModelID, effList=[10 ** lgeff], B0=1, kappa=0.5, compress=float(pase['compress']))
        elif len(ABC) == 2:
            """ Vary efficiency and B0"""
            (lgeff, lgB0) = ABC
            RModel = cbclass.RModel(RModelID, effList=[10 ** lgeff], B0=10 ** lgB0, kappa=0.5,
                                    compress=float(pase['compress']))
            print('#== Begin Processing task')
        elif len(ABC) == 3:
            """ Varies the standard model """
            (lgeff, lgB0, kappa) = ABC
            RModel = cbclass.RModel(RModelID, effList=[10 ** lgeff], B0=10 ** lgB0, kappa=kappa,
                                    compress=float(pase['compress']))
        elif len(ABC) == 4:
            """ Varies the standard model + detection probability """
            (lgeff, lgB0, kappa, survey_filter_pca_b) = ABC
            RModel = cbclass.RModel(RModelID, effList=[10 ** lgeff], B0=10 ** lgB0, kappa=kappa,
                                    compress=float(pase['compress']))
            surmodel = cbclass.SurModel(b=survey_filter_pca_b)
        elif len(ABC) == 6:
            (lgeff, lgB0, kappa, lgt0, lgt1, lgratio) = ABC
            RModel = cbclass.PreModel_Hoeft(RModelID, effList=[10 ** lgeff], B0=10 ** lgB0, kappa=kappa,
                                            compress=float(pase['compress']), t0=10**lgt0, t1=10**lgt1, ratio=10**lgratio)
            Rm = RModel
        elif len(ABC) == 7:
            (lgeff, lgB0, kappa, survey_filter_pca_b, lgratio, lgt0, lgt1) = ABC
            RModel = cbclass.PreModel_Hoeft(RModelID, effList=[10 ** lgeff], B0=10 ** lgB0, kappa=kappa,
                                            compress=float(pase['compress']), t0=10 ** lgt0, t1=10 ** lgt1,
                                            ratio=10 ** lgratio)
            Rm = RModel
            surmodel = cbclass.SurModel(b=survey_filter_pca_b)
        else:
            print('RunSurvey::main: model unknown')
            return


        """ Create survey """
        outfolder = '%s_%05i/' % (logfolder, RModelID)
        survey = cbclass.Survey(GClList, survey='%s' % (parfile.replace('.parset', '')),
                                emi_max=float(pase['RMSnoise']) * 1e-3 * 200,
                                cnt_levels=[float(pase['RMSnoise']) * 2 ** i for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]],
                                saveFITS=(ABC is None), savewodetect=suut.TestPar(pase['savewodetect']), dinfo=dinfo,
                                surshort='MUSIC2', Rmodel=RModel, outfolder=outfolder, logfolder=logfolder)
        survey.set_surmodel(surmodel)
        survey.set_seed_dropout()
    else:
        """ If you directly loaded a survey, just use its internal Rmodel """
        RModel = survey.Rmodel

    if verbose: print('Outfolder:', survey.outfolder)

    """=== Create a Task Cube, modify it & save it 
    The function of this Taskcube formerly was to keep track of all the computations that were already done in my multiproccecing version of ClusterBuster 
    Now it is outdated and isn't maintained anymore. After a BUG I didn't want to fix I decommisened this functionality
    """

    Taskcube = np.zeros((len(survey.GCls), len([RModel])))
    # A cube of all possible entries , efficiency is always fully computed and thus not in the Taskcube

    if suut.TestPar(pase['reCube']):
        Taskcube = np.load(logfolder + '/TaskCube.npy')
        smt.MergeSMT_simple(iom.unpickleObject(logfolder + '/smt'))
    if int(pase['reCL']) + int(pase['reRM']) > 0:
        for (GCl_ii, RModelID), status in np.ndenumerate(Taskcube[:, :]):
            if verbose: print('GCl_ii, RModelID:', GCl_ii, RModelID)
            if GClList[GCl_ii].mockobs.clid > int(pase['reCL']) and int(pase['reRM']):
                break
            else:
                Taskcube[GCl_ii, RModelID] = 1
    np.save(logfolder + '/TaskCube',
            Taskcube)  # also to be pickled/saved: Levels Cluster&TaskID --> B0, kappa, (z) -->  eff0
    """"""

    if verbose: print('#== Begin Processing task')

    """ This is the most important task! """
    while processTasks:
        processTasks, smt = DoRun((pase, survey), smt, verbose=verbose)

    print('RModelID %i of run %s finished' % (RModelID, surveyN))
    return survey

def main_ABC(params, parfile='MUSIC2_NVSS02_SSD.parset', Clfile='clusterCSV/MUSIC2', verbose=False):
    survey = main(parfile, ABC=params, Clfile=Clfile, verbose=verbose)

    """ MUSIC-2 """
    survey.dinfo = survey.GCls[0].dinfo
    survey.scatterkwargs = {"alpha": 0.15, "fmt": "^", "markersize": 7}
    survey.histkwargs = {"alpha": 0.15}
    survey = suut.AddFilesToSurvey(survey, savefolder=survey.outfolder, verbose=False, clusterwise=False)
    return survey


# from http://stackoverflow.com/questions/2553354/how-to-get-a-variable-name-as-a-string-in-python
def LoadSnap_multiprocessing(pase, realisations, Rmodel, getstring=False, verbose=False):
    smt = iom.SmartTiming()
    smt(task='LoadSnaps')   
           
    gcl = realisations[0]

    """ We load the radio cubes """

    strSn = (pase['snapfolder'] + 'SHOCKS_%05i/cluster.%05i.snap.%03i.shocks') % (gcl.mockobs.clid, gcl.mockobs.clid, gcl.mockobs.snap)

#    if suut.TestPar(pase['useMiniCube']):     # subset of snapshot, pickled
#        strSn   = strSn.replace('cluster.', 'clusterSUBSET.') 
#        if verbose: print('Loading snapshot:',strSn)
#        print('___ RunSurvey::LoadSnap_multiprocessing::',strSn,pase['useMiniCube'],suut.TestPar(pase['useMiniCube']))
#        snap    = iom.unpickleObject(strSn)
#        try:
#            snap    = iom.unpickleObject(strSn)
#        except:
#            with open('/data/ErrorLog.txt',"a") as f:
#                for gcl in realisations:
#                    f.write(strSn+'\n')   
#               
#    else: # original snapshot                          
#        if verbose: print('Loading snapshot:',strSn)
#        snap    = loadsnap.Loadsnap(strSn,headerc=pase['headerC'])    # original snapshot
        

    if  suut.TestPar(pase['useMiniCube']):   # asks: original or modified  snapshot?
         strSn   = strSn.replace('cluster.', 'clusterSUBSET.') 
    if verbose: print('Loading snapshot:',strSn)
    snap    = loadsnap.Loadsnap(strSn,headerc=pase['headerC'])                                    
                                    
                                    
    
    """psi and machfiles could become sharred (or global) arrays, but as they are quite small < 1MB, the impact on performance should be snmall"""
    PreSnap = radiomodel.PrepareRadioCube(snap, psiFile=pase['miscdata']+pase['PSItable'], machFile=pase['miscdata']+pase['DSAtable'])
    PreSnap = ( radiomodel.PiggyBagSnap(PreSnap[0]), PreSnap[1] )

    if getstring:
        return PreSnap, strSn, smt
    else:
        return PreSnap, smt




def varname(var):
  
  for k, v in list(locals().iteritems()):
         if v is var:
             a_as_str = k
  
  return a_as_str #[ k for k,v in locals().iteritems() if v is var][0]

  
def pool_wait(queues, limit, affix='', tmax = 1e3):
    tsleep = 0      
    waitT = 0.3

    while sum([q.qsize() for q in queues]) > limit and tsleep < tmax:  #stage 3 neclected
  
        if tsleep == 0:
             message = "[%s] Gonna sleep, because of " % (affix)
             stringlist = ['+%s(%i)'  % (varname(q),q.qsize()) for q in queues]
             string = ' '.join(stringlist)
             message += string
             message += " > %i" %(limit)
             print(message) 
        time.sleep(waitT) 
        tsleep += waitT
    print('pool_wait:', sum([q.qsize() for q in queues])) 
        
    if tsleep > 0:
        print("[%s] Slept for %.1f seconds. We don't want to shedule our memory to dead, do we?" % (affix, tsleep)   )
    
    return
    

def RadioAndMock_loaded(val, verbose=True):
    smt = iom.SmartTiming()
    smt(task='RadioAndMock_initialization')                
    (snapMF,  pase, realisations, survey) = val
    Rmodel = survey.Rmodel

    ##=== Stage II - DoMockObs
    if verbose:
        print('Start compiling MockObservations for further models of cluster #%5i snap #%3i with in total %i realisations.'
              % (realisations[0].mockobs.clid, realisations[0].mockobs.snap, len(realisations)))

    GClrealisations_return = []
    smt(task='Shed_DoMockObs_misc') 
    # This result of this computation is  independent of rotation and because of this was put here
    
    for kk, realisation in enumerate(realisations): # update information on the total radio power (at the rest frame frequency) in the simulational volume
        """ This is wrong and has to be fixed in the future!!!! 
        Currently, we make this code really SMELLY and hard to understand
        """
        (radiosnap, subsmt, poisson_factor) = radiomodel.CreateRadioCube(snapMF, Rmodel, realisation.mockobs.z_snap,
                                                         nuobs=pase['nu_obs'], logging=False)[0:3]
        smt.MergeSMT_simple(subsmt, silent=True)

        """ also possible: realisations[kk].Rvir       = radiocube[0].head['Rvir']"""
        if realisations[kk].M200.value == 0:
            try:
                realisations[kk].M200.value = radiosnap.head['M200']
            except:
                realisations[kk].Mvir.value = radiosnap.head['Mvir']
            realisations[kk].updateInformation(massproxis=True)

        """ Here we add the radio emission due to pre-existing electrons """
        if isinstance(Rmodel, cbclass.PreModel_Gelszinnis):
            randfactor = 10**np.random.normal(0, Rmodel.p_sigma, 1)
            realisation.PreNorm = randfactor*Rmodel.p0
            radiosnap.radiPre += realisations.PreNorm * radiosnap.radiPre
        elif isinstance(Rmodel, cbclass.PreModel_Hoeft):
            realisation.poisson_factor = poisson_factor

        """ Here we compute the volume weighted radio emission """   
        radiosum      = np.sum(radiosnap.radi)
        borders       = 2*realisations[0].R200*radiosnap.head ['hubble']/radiosnap.head ['aexpan']
        whereR200     = np.where(np.sqrt(np.power(radiosnap.pos[:,0], 2) + np.power(radiosnap.pos[:,1], 2)
                                         + np.power(radiosnap.pos[:,2], 2)) < borders/2)
        radiosum_R200 = np.sum(radiosnap.radi[whereR200])
        # update information on the total radio power (at the rest frame frequency) in the simulational volume
        realisations[kk].P_rest.value    = radiosum       # This is for a frequency differing from 1.4 GHz and an efficiency of 1, in PostProcessing.py we apply a further correction
        realisations[kk].Prest_vol.value = radiosum_R200  # This is for a frequency differing from 1.4 GHz and an efficiency of 1, in PostProcessing.py we apply a further correction

        if not suut.TestPar(pase['cutradio']):

            radiosnapUse = copy.deepcopy(radiosnap)
            if hasattr(radiosnapUse, 'radiPre'):
                radiosnapUse.radi += radiosnapUse.radiPre
                #print('Run_MockObs:: Ratio of PREs to total emission', (np.sum(radiosnapUse.radiPre)) / (np.sum(radiosnapUse.radi) + np.sum(radiosnapUse.radiPre)))

            radiocube = (radiosnapUse, Rmodel, survey)  # Rmodel is added to the tuple
        else:
            print('Beware, This is slow and should not be paralised!!!! This part is not implemented')
            radiocube = (radiosnap, Rmodel, survey)  # Rmodel is added to the tuple

        smt(task='Shed_DoMockObs')
        (nouse, subsmt, GClrealisation_used, Rmodel) = mockobs.Run_MockObs(radiocube, [realisation],
                                                                           saveFITS=survey.saveFITS,  savewodetect=survey.savewodetect,
                                                                           side_effects=True)
        GClrealisations_return += GClrealisation_used
        smt.MergeSMT_simple(subsmt, silent=True)
        
    return (GClrealisations_return, Rmodel), smt


def RadioAndMock(val, verbose=True):
    smt = iom.SmartTiming()
    smt(task='RadioAndMock_initialization')

    (pase, realisations, survey) = val
    Rmodel = survey.Rmodel

    PreSnap, smt_add = LoadSnap_multiprocessing(pase, realisations, Rmodel)

    if len(PreSnap) > 0:
        snapMF = PreSnap
        (radiosnap, subsmt) = radiomodel.CreateRadioCube(snapMF, Rmodel, realisations[0].mockobs.z_snap, nuobs=pase['nu_obs'])[0:2]
        smt.MergeSMT_simple(subsmt, silent=True) 

    """ This weird interresult comes from
        output.put( outp + (Rmodel,)) #Rmodel is added to the tuple
        (radiocube, subsmt, Rmodel) = stage1_out.get()
        stage1_list.append( ( radiocube, Rmodel, survey) )
    """
    radiocube = (radiosnap, Rmodel, survey) #Rmodel is added to the tuple


    ##=== Stage II - DoMockObs
    if verbose: print('Start compiling MockObservations for further models of cluster #%5i and snap #%3i with in total %i realisations.' % (realisations[0].mockobs.clid , realisations[0].mockobs.snap,  len(realisations)))


    smt(task='Shed_DoMockObs') 
    # This result of this computation is  independent of rotation and because of this was put here
    for kk, real in enumerate(realisations): # update information on the total radio power (at the rest frame frequency) in the simulational volume
             """ This is wrong and has to be fixed in the future!!!! """
             
             """ also possible: realisations[kk].Rvir       = radiocube[0].head['Rvir']"""
             if realisations[kk].M200.value == 0:
                 try:
                     realisations[kk].M200.value = radiocube[0].head['M200']
                 except:
                     realisations[kk].Mvir.value = radiocube[0].head['Mvir']
                 realisations[kk].updateInformation(massproxis=True)
             
    """ Here we compute the volume weighted radio emission """         
    radiosum  = np.sum(radiocube[0].radi)
    borders   = 2*realisations[0].R200*radiocube[0].head ['hubble']/radiocube[0].head ['aexpan']
    whereR200 = np.where( np.sqrt(np.power(radiocube[0].pos[:,0],2) + np.power(radiocube[0].pos[:,1],2) + np.power(radiocube[0].pos[:,2],2) ) < borders/2 )
    radiosum_R200 = np.sum(radiocube[0].radi[whereR200])

    # update information on the total radio power (at the rest frame frequency) in the simulational volume
    for kk,real in enumerate(realisations):
         realisations[kk].P_rest.value    = radiosum       # This is for a frequency differing from 1.4 GHz and an efficiency of 1, in PostProcessing.py we apply a further correction
         realisations[kk].Prest_vol.value = radiosum_R200  # This is for a frequency differing from 1.4 GHz and an efficiency of 1, in PostProcessing.py we apply a further correction

    smt(task='Shed_DoMockObs_misc')
    locations = [survey.outfolder]
    
    if not suut.TestPar(pase['cutradio']):
        radiocubeUse = radiocube
    else:  
        print('Beware, This is slow and should not be paralised!!!! This part is not implemented')
        radiocubeUse = None #radiomodel.PiggyBagSnap_cut(snap, radiocube[0], float(pase['cutradio'])),Rmodel,survey)]
    (nouse, subsmt, GClrealisations_used, Rmodel) = mockobs.Run_MockObs(radiocubeUse, realisations, saveFITS=survey.saveFITS, savewodetect=survey.savewodetect, writeClusters=True) #Mach=pase['Mach'], Dens=pase['Dens'],
    smt.MergeSMT_simple(subsmt, silent=True)

    return (GClrealisations_used,  Rmodel), smt


def RadioCuts(val, compradio=False):
    smt = iom.SmartTiming()
    smt(task='RadioAndMock_initialization')

    (pase, realisations, survey) = val

    Rmodel = survey.Rmodel
    snapMF, strSn, smt_add = LoadSnap_multiprocessing(pase,realisations,Rmodel, getstring=True)
    if compradio:   
        (snap,subsmt) = radiomodel.CreateRadioCube(snapMF, Rmodel, realisations[0].mockobs.z_snap, nuobs=pase['nu_obs'])[0:2]
        smt.MergeSMT_simple(subsmt, silent=True)
        
        smt(task='UpdateHeader')  
        realisations[0].z.value    = realisations[0].mockobs.z_snap
        realisations[0].Mvir.value = snap.head['Mvir']*snap.head['hubble']
        realisations[0].updateInformation()
        snap.head['M200'] = realisations[0].M200.value
                     
        """ Here we compute the volume weighted radio emission """         
        radiosum      = np.sum(snap.radi)    
        borders       = 2*realisations[0].R200*snap.head ['hubble']/snap.head ['aexpan']
        whereR200     = np.where( np.sqrt(np.power(snap.pos[:,0],2) + np.power(snap.pos[:,1],2) + np.power(snap.pos[:,2],2) ) < borders/2 )
        radiosum_R200 = np.sum(snap.radi[whereR200])

        # update information on the total radio power (at the rest frame frequency) in the simulational volume
        for kk,real in enumerate(realisations):
             realisations[kk].P_rest.value    = radiosum       # This is for a frequency differing from 1.4 GHz and an efficiency of 1, in PostProcessing.py we apply a further correction
             realisations[kk].Prest_vol.value = radiosum_R200  # This is for a frequency differing from 1.4 GHz and an efficiency of 1, in PostProcessing.py we apply a further correction
    else:
        (snap, MF) = snapMF

    smt(task='Shed_DoMockObs_misc')
    if suut.TestPar(pase['redSnap']): 
         # This is all part of reducing the data load. A shrinked version of the snapshot with less particles is saved in the custom format that is a fork of the gadget format
         redsnap = radiomodel.PiggyBagSnap_cut(snap, snap, float(pase['cutradio']), cutmask = MF )
         strSn   = strSn.replace('cluster.', 'clusterSUBSET.')
         strSn   = strSn.replace(pase['snapfolder'], pase['outf'])
         print(strSn)
         redsnap.savedata(strSn)    # original snapshot
    if compradio:
         print('np.sum(radiocube[0].radi)', np.sum(snap.radi), '-->', np.sum(redsnap.radi), 'i.e.', radiosum, radiosum_R200, 'pickled to', strSn)
         
    return (realisations, Rmodel), smt


def mupro_Output_NicePickleClusters( in_queue, output):
# Just pickles the relic for the highest efficiency, mechanism: It asks for the image   

    for well in in_queue:
        (Clrealisations,  Rmodel) = well
        outputMod = output #+ '_%05i/' % (Rmodel.id)
        iom.check_mkdir(outputMod + '/pickled')
        for GCl in Clrealisations:
           filename = 'MonsterPickle-Snap%02i' % (GCl.mockobs.snap)
           iom.pickleObject( (GCl, Rmodel), outputMod + '/pickled/', filename, append = True) 
    return


def DoRun(inputs, smt, verbose=False, countmax=500, countmax_relics=1500):
    """ Please mind that this procedure determines early if the number of detected relics becomes to large!"""  
    (pase, survey) = inputs

    count = 0
    count_relics = 0
    realisations = []
    realisations_list = [] # rotated realisations of one and the same cluster
    survey.GCls = sorted(survey.GCls, key= iom.Object_natural_keys)
    # This checks if the snapshot has to be loaded ...  I would use a pool of loaded, snapshots ... but this is cumbersome
    # Each snapshot is loaded only once ... If the number of models&rotations to be tested is high enough the produced overhead is low
    # This requires that the clusters are ordered due to cluster number and snapshot
    for gcl in survey.GCls:
        if (len(realisations) == 0 or gcl.mockobs.clid != realisations[0].mockobs.clid or gcl.mockobs.snap != realisations[0].mockobs.snap):
            if len(realisations) > 0:
                realisations_list.append(realisations)
            realisations = [gcl] 
        else: 
            realisations.append(gcl)   
            
    realisations_list.append(realisations)       

    for realisations in realisations_list:
            gcl = realisations[0]
            if verbose:
                print('Recognized ID %5i, snap %3i, #%i in cluster file' % (gcl.mockobs.clid, gcl.mockobs.snap, gcl.mockobs.id) )

            smt(task='LoadSnaps')   
            """ We load the radio cubes """
            strSn = (pase['snapfolder'] + 'SHOCKS_%05i/cluster.%05i.snap.%03i.shocks') % (gcl.mockobs.clid, gcl.mockobs.clid, gcl.mockobs.snap)
            """ SMELLY: This causes some issues, as there are two different 'load' versions for one and the same task """
            if suut.TestPar(pase['useMiniCube']):   # original snapshot
                strSn = strSn.replace('cluster.', 'clusterSUBSET.')
            if verbose:
                print('Loading snapshot:', strSn)
            snap = loadsnap.Loadsnap(strSn, headerc=pase['headerC'])

            smt(task='PrepareRadioCubes')   
            PreSnap = radiomodel.PrepareRadioCube(snap, psiFile=pase['miscdata']+pase['PSItable'],
                                                  machFile=pase['miscdata']+pase['DSAtable'])


            if 1==2:
                snap_phd = radiomodel.PrepareRadioCube(snap, psiFile=pase['miscdata'] + pase['PSItable'],
                                                      machFile=pase['miscdata'] + pase['DSAtable'], machmin=0.0)[0]

                print(type(snap_phd.mach))
                Mstat = np.asarray(snap_phd.mach)/1.045+1e-5
                Astat = np.asarray(snap_phd.hsml)/np.asarray(loadsnap.comH_to_phys(snap_phd.head))**2
                Rhostat = np.asarray(snap.rdow)
                Tstat = np.asarray(snap.udow)
                Psistat = np.asarray(snap.DSAPsi)

                bins = np.linspace(-1, 3, num=1300, endpoint=True)
                print(Mstat.shape, Tstat.shape, Psistat.shape)
                M_hist, bin_edges = np.histogram(np.log10(Mstat), bins=bins)
                A_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Astat, bins=bins)
                Rho_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Rhostat, bins=bins)
                T_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Tstat, bins=bins)
                Psi_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Psistat, bins=bins)
                Psi_A_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Psistat*Astat, bins=bins)
                Psi_A_T_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Psistat*Astat*Tstat**(1.5), bins=bins)
                Psi_A_T_Rho_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Psistat*Rhostat*Astat*Tstat**(1.5), bins=bins)
                Psi_A_Rho_hist, bin_edges = np.histogram(np.log10(Mstat), weights=Rhostat*Astat, bins=bins)

                alltogether = np.stack((M_hist, A_hist, Rho_hist, T_hist, Psi_hist, Psi_A_hist, Psi_A_T_hist, Psi_A_T_Rho_hist,Psi_A_Rho_hist), axis=1)
                print( len(snap.mach), alltogether.shape)

                for fr in realisations:
                    np.save("/data/TestForPhD/TestForPhD-%i" % (fr.ids), alltogether)

            PreSnap = (radiomodel.PiggyBagSnap(PreSnap[0], extended=True), PreSnap[1])
            if verbose:
                print(' ___ Particles in snap ___ :', PreSnap[0].radi.shape)  # DEBUGGING

            """ now we'll assign functions to the pool """
            inargs = (PreSnap, pase, realisations, survey)  #(radiocube, B0, kappa, z, strSn, interpolate)  # , int(TestPar(pase['default']))
            stage1_outSmall, smt_small = RadioAndMock_loaded(inargs, verbose=verbose)
            smt.MergeSMT_simple(smt_small)

            mupro_Output_NicePickleClusters([stage1_outSmall], survey.outfolder)
            """ pickles a bunch of single galaxy clusters, they are latter joined into a survey """
            """ It should also pickle the corresponding fits images (in case of a detection) """

          
            """ Block BEGIN: If you choose a poor model then you might end up with way more relic clusters then expected. 
                I presume that those both slows down the process and leads to working memory issues.
                To prevent this I stop the computation if the weighted relic count is way higher than the expected relic count.
                In this case I should let the procedure stop early and  give it the minimum score in the metric ...(also in the metric test)
            """
            for gcl_test in stage1_outSmall[0]:
                # print(survey.outfolder, gcl_test.stoch_drop(survey.seed_dropout))
                count_relics_now = len(gcl_test.filterRelics(**survey.relic_filter_kwargs))
                if count_relics_now > 0:
                    count += 1
                    count_relics += count_relics_now
                    #print(count, count_relics)
                if count > countmax or count_relics > countmax_relics:
                    print('Because the number of relics already now is larger than %i the function DoRun() is terminated.' % countmax_relics)
                    return False, smt
            """  Block END """

            realisations = []  # Clear rotations of clusters, so that append, can work for a new try, ...
            # formerly: del realisations


    if verbose: print("#=== Proccesses completed. Please check for the data output.'")
    if verbose: smt(forced=True)

    return False, smt


def create_ClusterLibrary(snapfolder='/radioarchive/MergerTrees/Clusters/snaps/', Clfile='clusterCSV/MUSIC2', copyfolder='/data2/MUSIC-2/',
                          clusters =range(1,283), snaps =range(3,18), add=''):
    """ This is an aunciiliary function that should be run once to create a file that lets you map R200 and M200 to each snapshot.
        This is done, because in the current implementation of the ClusterBuster_pythonClass.py the galaxycluster computing these quantities is quite costly
        The cost is that the taken R200 and M200 differ from the real values at zobs != zsnap,
        
        BUG and SMELL: THIS routine however is not used anymore (?), as the Radiocuts computes its own M200 values, before it compresseses the snapshots!"""
    
    from shutil import copyfile
    import os.path

                
    """ UTILITY FOR THESIS           
    for snapID in [0,1,2,3,4,5,6]:            
        for clid in clusters:
                folderSn = '%s%s/SHOCKS_%05i/' % (snapfolder, '', clid)
            
                strSn    = '%s%s/SHOCKS_%05i/cluster.%05i.snap.%03i.shocks' % (snapfolder, '', clid, clid, snapID)
                if os.path.isfile(strSn):
                   # print('rm %s'  % (strSn))
                    os.system('rm %s'  % (strSn))           
    return 0
    UTILITY FOR THESIS END """     

    MissList     = []  
    index        = 0
    numberOfRows = 0
    for snapID in snaps:
        Outlist = []
        for clid in clusters:
            for adds in [add]:
                strSn = '%s/SHOCKS_%05i/cluster.%05i.snap.%03i.shocks' % (snapfolder, clid, clid, snapID)
                if os.path.isfile(strSn):
                    numberOfRows += 1
                    
    # create dataframe
    df = pd.DataFrame(index=np.arange(0, numberOfRows), columns=('clID', 'snapID', 'redshift', 'M200', 'originalpath') )


    for snapID in snaps:
        Outlist = []
        for clid in clusters:
            for adds in [add]:
                """ formerly this was ['','Additional2'] because different snapshots were stored in different folders """
                strSn = '%s/SHOCKS_%05i/cluster.%05i.snap.%03i.shocks' % (snapfolder, clid, clid, snapID)
                print (strSn)
                if os.path.isfile(strSn):
                    snap  = loadsnap.Loadsnap(strSn)    
                    z     = 1/snap.head['aexpan']-1
                    cl = cbclass.Galaxycluster_simulation("",  0, z=z, Mvir=snap.head['Mvir']*snap.head['hubble']) 
                    # ID SNAP z M200  
                    df.loc[index] = [clid, snapID, z, cl.M200.value, strSn]    
                    index += 1
                    print(strSn + ' was loaded')
                    
                    
                    if copyfolder is not None:
                        folderSn = '%s/SHOCKS_%05i/' % (snapfolder, clid)
                        folderDe = '%s/SHOCKS_%05i/' % (copyfolder, clid)
                        strSn = folderSn + 'clusterSUBSET.%05i.snap.%03i.shocks.pickle' % (clid, snapID)
                        strDe = folderDe + 'clusterSUBSET.%05i.snap.%03i.shocks.pickle' % (clid, snapID)
                        print(strSn)
                        if os.path.isfile(strSn):
                            print(strSn + ' was recognized and is copied to %s ' % (copyfolder))
                            iom.check_mkdir('%s%s/SHOCKS_%05i/' % (copyfolder, adds, clid))
                            copyfile(strSn,strDe)  
                            
                    if index % 20 == 0:
                        df.to_csv('%s_allclusters.csv' % (Clfile))     
                else:
                    print('  skipping misslist step')
                    MissList.append('%5i %4i' % (clid, snapID))

        if index > 0:
            df.to_csv('%s_allclusters.csv' % (Clfile))
#            thefile = open('%s_z%03i.dat' % (Clfile, z*100), 'w') 
#            for item in Outlist:
#                thefile.write("%s\n" % item)
#            thefile.close()
            
    if len(MissList) > 0:
        thefile = open('%s-missing.dat' % (Clfile), 'w') 
        for item in MissList:
            thefile.write("%s\n" % item)
        thefile.close()   

    return 0
    

def ReloadSurvey(survey=None,parfrom='MUSIC2COOL_NVSS_SSD.parset', parto='MUSIC2_NVSS02_SSD.parset', index=2, index_new=3,  reform=True):
    """ Reloads and runs a single survey and the selected clusters again - after running it saves it into another directory
        This is usefull to test other survey configurations
        
        reform: Allows to reform the snapshot parameter format of the galaxy cluster. This is done to allow a direct comparison between 
         the runs with and without cooling
    """
    if survey is None:
        surveyname = parfrom.replace('.parset', '') + '_%05i' % (index)
        savefolder = '/data/ClusterBuster-Output/'
        print(parfrom, savefolder, surveyname)
        survey = iom.unpickleObject('%s/%s/pickled/Survey' % (savefolder, surveyname))
        
    """ reform is [-1,0,1] and is used to convert from the different snapshot numbering systems"""
    survey.name      = parto.replace('.parset', '') + '_%05i' % (index_new)
    survey.outfolder = savefolder + survey.name
    survey.logfolder = savefolder + survey.name

    
#    pase, _, _, _, _ = suut.interpret_parset(parfile, repository=workdir+'/parsets/', relative=(ABC is not None))
#    misslist     = np.loadtxt('/home/jakobg/lib/ClusterBuster/surveysim/music2/AllMUSIC-AGN-missing.dat')   #pase['miscdata']+pase['missFile'])
    misslist     = np.loadtxt('/home/jakobg/lib/ClusterBuster/surveysim/music2/clusterCSV/MUSIC2-missing.dat')   #pase['miscdata']+pase['missFile'])
    fitleredGCls = []
    missing      = False
    """Thesis plot, DEVLOPMENT"""
    """Thesis plot"""
    
    for GCl in survey.GCls:
        GCl.relics = []  #do not consider any previous detections
        if reform: 
            """ Changes that format from compatible with the radiative simulation to the adiabatic simulation """
            GCl.mockobs.snap -= 3     #or +3
            GCl.mockobs.headerc = 'Mpc' #or 'kpc'
            GCl.mockobs.snapfolder = '/data2/MUSIC-2/'
            GCl.mockobs.xrayfolder = '/radioarchive/MergerTrees/Clusters/'

        ids  = GCl.mockobs.clid
        snap = GCl.mockobs.snap
        """ Exclude missing snapshots """
        for miss in misslist:
            if (miss[0] == ids) and (miss[1] == snap):
                missing=True

        if not missing:
            fitleredGCls.append(GCl)
        else:    
            missing=False
    print(len(sorted([gcl.mockobs.id for gcl in survey.GCls])))
    survey.GCls = fitleredGCls
    
    if reform:
        parfile = parto 
    else:
        parfile = parfrom

    """ DEVELOPMENT  HOEFT_model
    survey.Rmodel = cbclass.PreModel_Hoeft(0, effList=[3e-6], B0=2.0, kappa=0.4)  #
    survey.Rmodel.t0      = 0.055   # Minimal time since reacceleration
    survey.Rmodel.t1      = 0.55    # Maximal time since reacceleration
    survey.Rmodel.n0      = 1e-6   # Number density of accretion shocks
    survey.Rmodel.n1      = 1e-2   # Number density of 'core'
    survey.Rmodel.ratio   = 0.05   # Initial normalisation PRE and thermal at shockfront
    """

    survey = main(parfile, workdir=None, ABC=None, verbose=False, survey=survey)
    time.sleep(16)
    suut.AddFilesToSurvey(survey, savefolder=survey.outfolder, verbose=False,clusterwise=True) 

    return 0

def copy_ClusterOuts(snapfolder = '/data/ClusterBuster-Output/', copyfolder = '/data2/MUSIC-2/'):
    
    import os.path
    import shutil

    runs =range(100000) # list of ...
    copyfolder = snapfolder
    for run in runs:
            for adds in ['']:
                """ formerly this was ['','Additional2'] because different snapshots were stored in different folders """
                folderSn = '%s/MUSIC2_NVSS02_SSD_%05i' % (snapfolder, run)
                folderDe = '%s/abcpmc_MUSIC_NVSS02_Run_02/MUSIC2_NVSS02_SSD_%05i' % (copyfolder, run) 
                os.system('mv %s.txt %s.txt' % (folderSn, folderDe))
                try:
                   shutil.copy ('%s.txt' % folderSn, ' %s.txt' % folderDe)
                except:
                   print('%s.txt' % folderSn, ' %s.txt' % folderDe, 'could not be copied')
    return 0


def test_main_ABC():
    #test = main_ABC([-7, -1.0, -2.0, -1.0, 0.0, -4.0], verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset
    test = main_ABC([-7, 0.0, 0.5], parfile='MUSIC2_NVSS01_SSD.parset', Clfile='clusterCSV/MUSIC2', verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset
    print("main_ABC returns:", type(test), test)

if __name__ == "__main__":
#    create_ClusterLibrary(snaps=range(0,15), copyfolder=None)
#    create_ClusterLibrary(snaps=range(3,18),snapfolder='/radioarchive/MergerTrees/WithCoolSfrAgn/snaps/', Clfile='clusterCSV/MUSIC2-AGN', copyfolder=None)

#    Create smaller snapshots
#    copy_ClusterOuts()
#    main('MUSIC2_prep.parset'  , Clfile = 'clusterCSV/AllMUSIC'    , verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset  
#    main('MUSIC2COOL_prep.parset', Clfile = 'clusterCSV/AllMUSIC-AGN', verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset        

#    create_ClusterLibrary(snapfolder = '/radioarchive/MergerTrees/WithCoolSfrAgn/snaps/', copyfolder = '/data2/MUSIC-2-AGN/')
#    create_ClusterLibrary(snapfolder = '/radioarchive/MergerTrees/Clusters/snaps/', copyfolder = '/data2/MUSIC-2-correctedM/')


# ABC run
#    main('%s.parset' % (args.par), Clfile = 'clusterCSV/AllMUSIC-AGN', verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset   

# DEBUGGING

#    snap   = loadsnap.Loadsnap('/data2/MUSIC-2/snaps/SHOCKS_00001/clusterSUBSET.00001.snap.004.shocks')      #,headerc='Mpc'
#    snap   = loadsnap.Loadsnap('/radioarchive/MergerTrees/Clusters/snaps/SHOCKS_00001/cluster.00001.snap.004.shocks')      #,headerc='Mpc'
#    print( snap.head)
#    print( np.mean(snap.pos[:,0]), np.mean(snap.pos[:,1]), np.mean(snap.pos[:,2]), snap.head )


#    main_ABC([-5,0,0]) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset
#        index
   
#    create_ClusterLibrary(snapfolder = '/data/Threehundrets/', Clfile = 'clusterCSV/Threehundrets', copyfolder = None, clusters =range(1,2), snaps  =range(86,129)) 
#    survey = main('Threehundrets_prep.parset', Clfile = 'clusterCSV/Threehundrets', verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset 
#    survey = main('Threehundrets.parset', Clfile = 'clusterCSV/Threehundrets', verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset
#    suut.AddFilesToSurvey(survey, savefolder=survey.outfolder, verbose=False,clusterwise=True)

#    survey = main('ShockTest_prep.parset', Clfile = 'clusterCSV/ShockTest', verbose=True) #FullRun_testNuza3.parset') #NVSS_Berlin00C.parset
    test_main_ABC()