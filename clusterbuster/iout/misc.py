#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 2017

@author: jakobg
"""

from __future__ import division,print_function

import os
import numpy as np
import time



#============== Smart timing to subject classes


class SmartTiming():
  """
  Class to allow timing of specific task blocks within the code. with 
  output about the total add taskspecific timecosnumption (real time, and cpu usage).
  Intermediate complexity.
  """

  def __init__(self, rate=5e4, logf=False): #, ini=False
    """
    decription
    """
    t                =  time.time() 
    tcpu             =  time.clock() 
    self.start       = [t,tcpu]
    self.now         = [t,tcpu]
    self.lastadd     = [t,tcpu]
    self.lastout     = t
    self.tasks       = [ ['unspecified',0,0] ]  #[Task name, processtime, cputime]
    self.tcurrent    = 'unspecified'
    self.rate        = rate
    self.logf        = logf
    self.ini         = True   
    self.pos         = 0
    
    
  def update_pos(self, pos):
    self.pos = max(self.pos,pos)
    
    
  def output(self, forced=False):
  
     if (self.now[0] - self.lastout > self.rate) or forced:
       print( '#==== Timing report BEGIN \n With a total computing time of %.1f seconds the following tasks took so much of the total time:' %(self.now[0]-self.start[0]) )
       for task in self.tasks:
          print( "Task '%24s' ran %5.1e seconds accounting for %7.3f %% of the total computing time. Average CPU usage was %7.3f %%." % ( task[0], task[1], 1e2*task[1]/(self.now[0]-self.start[0]),1e2*task[2]/task[1] ))
        
       print( '#==== Timing report END \n ')
       self.lastout = self.now[0]
       if self.logf:

         if self.ini:
           with open(self.logf, 'w') as file:
            line = '%17s ' % ('Ncl     t_total [s]')
            for task in self.tasks:
               line += '%24s  ' % (task[0])
            file.write(line+'\n')
            self.ini = False

         with open(self.logf, 'a') as file:
           line = '%7i %10.1f  ' % (self.pos, self.now[0]-self.start[0])
           for task in self.tasks:
             line += ' %5.1e %7.3f %7.3f' % (task[1], 1e2*task[1]/(self.now[0]-self.start[0]), 1e2*task[2]/task[1])
           file.write(line+'\n')
    
  def __call__(self, task='unspecified', forced=False, pos=0):
    self.now     = [time.time(),time.clock()]
    tdiff        = [x - y for x,y in zip(self.now,self.lastadd)]   
    try: 
      self.tasks[ [x[0] for x in self.tasks].index(self.tcurrent)][1] += tdiff[0] 
      self.tasks[ [x[0] for x in self.tasks].index(self.tcurrent)][2] += tdiff[1]
    except:
      self.tasks.append( [self.tcurrent, tdiff[0], tdiff[1]] )
      #print( self.tasks[self.tasks[0][:].index(self.tcurrent)]
       
    self.update_pos(pos)
    self.output(forced=forced)
    self.lastadd = self.now
    self.tcurrent = task
    
  def MergeSMT(self, smt):   #Merges to smt tasks ... very bulky however

    for jj,task in enumerate(smt.tasks): 
        	try: 
        	  self.tasks[ [x[0] for x in self.tasks].index(task)][1] += task[jj][1]
        	  self.tasks[ [x[0] for x in self.tasks].index(task)][2] += task[jj][2]
        	except:
        	  self.tasks.append( task )  


  def MergeSMT_simple(self, smt, silent = True):   #This implies that the smt tasks are complete! It wont add timings of new tasks
  
    #print( 'Check:', [t[0] for t in self.tasks],[t[0] for t in smt.tasks]
    for jj,task_b in enumerate(smt.tasks):
        found = False
        for ii,task in enumerate(self.tasks):  
            if task[0] == task_b[0]: # Compare task names
                # print( task[0], self.tasks[ii][1], smt.tasks[jj][1], self.tasks[ii][2], smt.tasks[jj][2]
    
                self.tasks[ii][1] += smt.tasks[jj][1] # Add up CPU time
                self.tasks[ii][2] += smt.tasks[jj][2] # Add up System Time
                found = True
                continue     # Exit loop
	if not found:
	   self.tasks.append( task_b )  
	   if not silent: print( 'Task not known!!! Task list is getting appended', task_b)
               

    
  
 
import cPickle as pickle
#import dill as pickle   # To make lambda function pickleable, use with caution
# See http://stackoverflow.com/questions/4529815/saving-an-object-data-persistence-in-python for an even more sophisticated alternative

def pickleObject(obj, location, oname, append = False):

  check_mkdir(location)  
  commands = ['wb','ab']
  with open(location+oname+'.pickle', commands[int(append)]  ) as handle:
#      print('misc::pickleObject()::', location+oname+'.pickle') #DEBUGGING
      pickle.dump(obj, handle, -1)
      
  return None    
    
def pickleObject_old(obj, location, append = False):
    ''' With the check to create the directory, but with the full path '''
    check_mkdir(os.path.dirname(location)) ## directory of file

   
    commands = ['wb','ab']
    with open(location+'.pickle', commands[int(append)]  ) as handle:
        pickle.dump(obj, handle, -1)
      
    return None    
  
def unpickleObject(location):
  
  with open(location+'.pickle', 'rb') as handle:
      obj = pickle.load(handle)
  return obj
    
    
# from Lutz Prechelt: http://stackoverflow.com/questions/20716812/saving-and-loading-multiple-objects-in-python-pickle-file  
def unpickleObjectS(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break
    
def load_file_clusters ( filename ) :
    with open( filename, 'r' ) as f:
        line     = f.readline()
        colnameL = line.lstrip('#').rstrip('\n').split()
        line     = f.readline()
        line     = f.readline()
        dtypeL   = line.lstrip('#').rstrip('\n').split()
    f.close()

    NCol     = len(colnameL)
    TabShape = zip(colnameL,dtypeL)

    valueLL = []
    for line in open( filename, 'r' ) :
        tokenL = line.strip('\n').split()
        if len(tokenL) < NCol  or  tokenL[0][0]=='#' :
            continue
        valueLL.append(tuple(tokenL))

    table = np.array( valueLL, dtype=TabShape )

    return table


LabelD = {
           'flux':'flux [mJy]',
     'a_rel_iner':'angel relative relic orientation to position  [deg]',
       'iner_rat': 'ratio of moments of inertia'
}
    

def read_para_list( ParaFile, log=False ) :
    
   if log:
       print( ' reading parameters')
       print( '    filename: %s' % ParaFile)

   scL = []

   for nl, line in enumerate(open( ParaFile )) :
     
      token = line.split('|')
      
      if len(token) < 3 or token[0][0]=='#' :
         continue

      for i, el in enumerate(token) : 
         token[i] = el.rstrip('\n').rstrip(' ').lstrip(' ')
         
      if nl == 0:
         keyL  = []
         typeL = []
         for el in token[0:] :
            tok2 = el.split()
            typeL.append(tok2[0])
            keyL.append( tok2[1])

         valueStore = []
         for key in keyL :
            valueStore.append('')
      else:
         dict = {}
         for i, el in enumerate(token[0:]) :
            key = keyL[i]
            value = el   #without value store
            if   typeL[i] == 'i' :
               dict[key] = int(value)
            elif typeL[i] == 'f' :
	       try:
                 dict[key] = float(value)
               except:  
                 dict[key] = 0
            elif typeL[i] == 'b' :
               dict[key] = (value == 'True')
            else :
               dict[key] = value

         scL.append(dict)

   if log:
       for i, el in enumerate(scL) :
          print( '   %3i' %i, el )
       print( ' reading of parameters done' )
             
   return scL
   
   
def parset2dict(parfile, relative=False):
    
  if relative:
      import os
      parfile = os.getcwd() + parfile
  with open(parfile) as f:
      lines = f.readlines()
#  import re 
  new_dict = dict()
  for ii, line in enumerate(lines[:]):
    if line[0] not in ['#'] and line.replace(" ","").replace('\n',''):   
      line  = line[:-1]
      parts = line.split('#')
      args  = parts[0].split(' ')  # Splits arguments by ' '
      args  = filter(None, args)   # Removes empty string args
      #print( args 
      new_dict[str(args[0])] = args[2]   # We should have a list of the args now [parameter,=,value] ... comments not included
        
  return new_dict
   

#import re 
  #====

def CoordinateToFloats(COO):
    sign = 2*int(COO[3] >= 0)-1  
    RA   = 15*(COO[0]  +       COO[1]/60  +      COO[2]/3600)
    Dec  = 1 *(COO[3]  +  sign*COO[4]/60  + sign*COO[5]/3600)
    return (RA, Dec)

def CoordinateToPixel(COO, spixel, center, pixref, pixels = True):
  

     COO  = COO.astype(float)
     RA, Dec = CoordinateToFloats(COO)
     
     if pixels:
         PixelCnt = np.zeros( (2), dtype=float)
         PixelCnt[0] =  (np.cos(Dec*np.pi/180)*(RA  - center[0]) / spixel[0]) + pixref[0]   # RA
         PixelCnt[1] =  (                      (Dec - center[1]) / spixel[1]) + pixref[1]   # Dec
         return PixelCnt
     else:
         return np.array( (RA, Dec) )
  

def ContourToList(COOcnt, spixel, center, pixref, pixels=True):
  COOcnt = (COOcnt.replace(':',',')).split(',')
  COOcnt = np.asarray(COOcnt)

  COOcnt  = np.reshape(COOcnt, (COOcnt.shape[0]/6, 6))
  
  if pixels:
     PixelCnt = np.zeros( (COOcnt.shape[0], 2), dtype=int)
  else:
     PixelCnt = np.zeros( (COOcnt.shape[0], 2), dtype=float)
     
  for ii,COO in enumerate(COOcnt):
     pixel          =  CoordinateToPixel(COO, spixel, center, pixref, pixels) 
     if pixels:
         PixelCnt[ii,0] =  int (pixel[0]+0.5)              
         PixelCnt[ii,1] =  int (pixel[1]+0.5)   
     else:
         PixelCnt[ii,0] =       pixel[0]*1.0           
         PixelCnt[ii,1] =       pixel[1]*1.0   
  return np.expand_dims(PixelCnt, axis=1) 

  
  
# from   astLib import astCoords      would be quite handy here


#from astropy import units as u
#from astropy.coordinates import SkyCoord
  
def J2000ToCoordinate(J2000):
  
  J2000 = J2000.replace("J2000 ","")
  
  # astropy solution
  #c = SkyCoord('00h42m30s', '+41d12m00s', frame='icrs')
  #c = SkyCoord(J2000[0], J2000[1], frame='icrs')  
  #print( c
  
  try: 
    COO = float(J2000)
  except:
    replace = ['h', 'm', 's', 'd']
    for r in replace:
       J2000 = J2000.replace(r,' ')
    COO = str(J2000)

    COO = np.asarray(filter(None, COO.split(' '))) 
  return COO

  
  
def readDS9regions(regfile, spixel, center, pixref, pixelcoords=True):

   contours     = []
   contourinfos = []
   with open(regfile) as line:
     for l in line:
       if 'polygon' in l:

         textpos     =   l.find('text={') 
         if textpos > -1: 
           contourinfo =   (l[textpos+6:].replace('}\n','')).replace(" dash=1","") 
         else:
           contourinfo    = ('')

         l       = l.split(' ')
         l_clear = (((l[0].replace("polygon(", "")).replace("\n","")).replace(")","")).replace("dash=1","")
         contours.append( ContourToList(l_clear, spixel, center, pixref, pixels = pixelcoords) )  
         contourinfos.append(contourinfo)
   return (contours, contourinfos) #Contour list of numpy arrays(Npoints, 1, 2)
   
   
   
def plot_smt(folder, smtlog,log=False, rel=False):  ##plot_smt('../Analysis-MUSIC-2/Output-JG-N40000-eff0.2/smt.log' )
  import re
  
  with open(folder + smtlog, 'r') as f:
    first_line = f.readline()
    header = re.split(r'[\s]\s*',first_line)
    
   
  Ntasks = (len(header)-3)

  data = np.loadtxt(folder + smtlog,skiprows=1) #In fact skiprows=1, but currently (some of the) first rows could be bugged due to missing entries
  print( Ntasks, data.shape )
  
  import matplotlib.pyplot as plt
  
 #(task[1], 1e2*task[1]/(self.now[0]-self.start[0]), 1e2*task[2]/task[1])
 
  xunit = 1000
  print( data )
  x = data[:,0]/xunit  #in k realisations
  style =  ['-',':']
  for ii in np.arange(0,Ntasks-1):

    y = data[:,ii*3+2]
    ylabel = '\\mathrm{cpu\\,time \\,[s]}' 
    linewidth = 1 if ii<8 else 3
    alpha     = 1 if ii<8 else 0.3
    
    if rel: 
       y= y/(x*xunit)
       ylabel += '\\mathrm{\\,per\\,realisation}'
    if log: 
       y = y+1e-4
       plt.yscale('log')

    linestyle =   style['sub'  in header[ii+3]]
    plt.plot(x, y,label=header[ii+3],linewidth=linewidth,linestyle=linestyle,alpha=alpha)
    #plt.scatter(x, y,label=header[ii+3], alpha=alpha, s=1)
    
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=0,
           ncol=4, mode="expand", borderaxespad=0., fontsize=5)

  if log:      
       lims = plt.ylim()
       plt.ylim([-2, lims[1]]) 
       ylabel = '$\\mathrm{log_{10}}(' + ylabel + ')$'
           
  #plt.show()
  lims = plt.xlim()
  plt.xlim([0, lims[1]]) 
  plt.xlabel('# of realisation in units of thousands')
  plt.ylabel(ylabel)
  plt.savefig(folder + 'smt.pdf')
  

  
def str2list(string):
  
  try:
    string = string[1:-1]
    lis    = string.split(',')
  except:
    print( 'String:',string,'could not be parset into a list' )
  
  return lis
  
def createlist(Prange,Nbin,interpol='lin'):
  
  if 'lin' in interpol:
    lis = np.linspace(float(Prange[0]), float(Prange[1]), int(Nbin), endpoint=True)  #
  elif 'log' in interpol:
    lis = np.logspace(np.log10(float(Prange[0])), np.log10(float(Prange[1])), int(Nbin), endpoint=True)  #
  else:
    print( "createlist(Prange,Nbin,interpol='lin'): ????" )
    
    lis = []
    
  return lis
  
#==== Natural number sorting, taken from http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
    
def Object_natural_keys(o): #, reverse = False
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    #if reverse: ll= ll[::-1]check_mkdir
    return [ atoi(c) for c in re.split('(\d+)', o.name) ]
    
    
def memory_usage_psutil(process): #taken from http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
    # return the memory usage in MB
    #import psutil
    #process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20) #process.get_memory_info()[0] / float(2 ** 20)
    #print( '!!!!!process.memory_info().rss', process.memory_info().rss/ float(2 ** 20) #seems to be same
    return mem  
    
    
def check_mkdir(filename): # taken from  http://www.stealthcopter.com/blog/2009/09/python-making-multi-depth-directories/

#    # Variant A
#	folder=os.path.dirname(filename)
#	if not os.path.exists(folder):
#		os.makedirs(folder)
    #Variant B
    try:
        os.stat( filename )
    except:
        os.makedirs( filename ) 
    
  
