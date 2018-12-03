#!/usr/bin/env python

# Modified by Jakob Gelszinnis

"""

 What is the standard of this data format anyhow. If we have a load function, do we also have save function other than .pickle 
 ... the issue of pickling is that it leads to very strong coupling.

"""



from __future__ import division,print_function

import os
import errno
import struct
import sys
import numpy as np
import clusterbuster.constants as cgsconst

const = cgsconst.ConstantsCGS()

def get_identi (f, end, verbose = False ) :
   bytes  = f.read( 4 )
   if len(bytes) != 4 :
      if verbose: print(' get_identi: no further valid identifier found')
      return None
   identy = ''.join( struct.unpack( end+'4c', bytes ) ) 
   return identy 
   
   
def Loadsnap(strSn, transform=True, headerc='Mpc'):
  
  #====  Hoeft_radio/q_mach_machr_table.txt for spectral index 
    snap = Snapshot(strSn, transform=transform, headerc=headerc)  #, radio_name = strRa
    snap.loaddata()
    snap.head['gamma'] = 5./3.  # add adiabatic expansion factor

    #==== psiFile   for psi factor; machfile for mach-numbers conversion factors
    return snap

def savesnap(strSn, savefolder, transform=True, headerc='Mpc'):
  
    snap = Snapshot(strSn, transform=transform, headerc=headerc)  #, radio_name = strRa
    snap.savedata(savefolder) 


def load_snap(snap, verbose=False):
    
   if verbose:
     print(' load snapshot <%s> :' % snap.name)
   
   try :
      f = open( snap.name, 'rb' )
   except IOError:
      print('Fatal Error: apparently %s does not exist' % (snap.name))
      print('... I better stop now !!!')
      sys.exit()
      
   end = '<' # args.endian
      
   # check FileInit
   FileInit = struct.unpack( end+'II',  f.read( 2*4 ) )
   if FileInit != (1,2) :
      print('Fatal Error: FileInit != (1,2)')
      print('(possibly something wrong with endian)')
      print('... I better stop now !!!')
      sys.exit()
      
   # read header
   NumPart, Mvir, Rvir, Xc, Yc, Zc, aexpan, hubble = struct.unpack( end+'lfffffff', f.read( 8*4 ) )
   """DEVELOPMENT, allows to use data that were already transformed"""
   if snap.headerc == 'Mpc':
       Xc *= 1e3   #  Mpc -> kpc
       Yc *= 1e3   #  Mpc -> kpc
       Zc *= 1e3   #  Mpc -> kpc

   elif snap.headerc != 'kpc':
       print('Snapshot header of central coordinates: UNIT "%s" unkown!' % (snap.headerc))
       
   if verbose:
      print('   NumPart = %6i'         % NumPart )
      print('   Mvir    = %.3e Msun/h' % Mvir    )
      print('   Rvir    = %.3e kpc/h'  % Rvir    )
      print('   Xc      = %.3e kpc/h'  % Xc      )
      print('   Yc      = %.3e kpc/h'  % Yc      )
      print('   Zc      = %.3e kpc/h'  % Zc      )
      print('   axpan   = %.3f    '    % aexpan  )
      print('   hubble  = %.3f    '    % hubble  )
   f.read( 512 - 10*4 )
   
   snap.head = { 'hubble':hubble, 'aexpan':aexpan, 'Xc':Xc, 'Yc':Yc, 'Zc':Zc, 'NumPart':NumPart, 'Mvir':Mvir, 'Rvir':Rvir }
   
   while True:
   
      identi = get_identi ( f, end )
      if identi == None :
         break
      
      if verbose: print(' following identifer found: <%s>' % identi)
   
#      print('_',identi)
      if identi == 'POS ':
         snap.pos  =  np.fromfile(f, dtype=end+'f', count=3*NumPart ).reshape( (NumPart, 3 ) )
         if snap.transform:
             snap.pos[:,0] -= Xc
             snap.pos[:,1] -= Yc
             snap.pos[:,2] -= Zc 
      elif identi == 'VEL ' :
         snap.vel  =  np.fromfile( f, dtype=end+'f', count=3*NumPart ).reshape( (NumPart, 3 ) )
      elif identi == 'MGAS' :
         snap.mgas =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'ID  ' :
         snap.id   =  np.fromfile( f, dtype=end+'I', count=NumPart )
      elif identi == 'HSML' :
         snap.hsml =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'U   ' :
         snap.u    =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'RHO ' :
         snap.rho  =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'ENDT' :
         snap.endt =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'NORM' :
         snap.norm =  np.fromfile( f, dtype=end+'f', count=3*NumPart ).reshape( (NumPart, 3 ) )
      elif identi == 'UUP ' :
         snap.uup  =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'UDOW' :
         snap.udow =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'RUP ' :
         snap.rup  =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'RDOW' :
         snap.rdow =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'DIVV' :
         snap.divv =  np.fromfile( f, dtype=end+'f', count=NumPart )
      elif identi == 'MACH' :
         snap.mach =  np.fromfile( f, dtype=end+'f', count=NumPart )
      else :
         print('Identifier is not recognized'  ) 
         break
      
   f.close()
   if verbose:
      print(' load snapshot ... done')

   return 



def save_snap ( snap, filename, verbose=False) :
    
   if verbose:
       print(' save snapshot <%s> :' % snap.name)
     
#   try :

   if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
   fwrite = open( filename, 'w' )
   end = '<'

   fwrite.write(struct.pack( end+'II',  1,2 ))
      
   """ The snapshot is altered, i.e. this is smelly codeas it should just save a snapshot. In our case only not of relevance because, we either save the snapshots or create radio maps."""
   if snap.transform:
         snap.pos[:,0] += snap.head['Xc']
         snap.pos[:,1] += snap.head['Yc']
         snap.pos[:,2] += snap.head['Zc']
    
   if snap.headerc == 'Mpc':
       snap.head['Xc'] /= 1e3   #  Mpc -> kpc
       snap.head['Yc'] /= 1e3   #  Mpc -> kpc
       snap.head['Zc'] /= 1e3   #  Mpc -> kpc
   elif snap.headerc != 'kpc':
       print('Snapshot header of central coordinates: UNIT "%s" unkown!' % (snap.headerc))
   
 

   snap.head['NumPart'] = snap.hsml.shape[0]
   fwrite.write(struct.pack( end+'lfffffff', snap.head['NumPart'], snap.head['Mvir'], snap.head['Rvir'], snap.head['Xc'], snap.head['Yc'], snap.head['Zc'], snap.head['aexpan'], snap.head['hubble'] ) )

#   f.read( 512 - 10*4 )
   fwrite.write( '{0:0472b}'.format(0) )
   
   pairs  = [ ['POS ',lambda x:  x.pos],
              ['VEL ',lambda x:  x.vel],
              ['NORM',lambda x:  x.norm] ]
   for ident,lambd in pairs:
#       print(struct.pack( end+'s', ident), struct.pack( end+'c', ident[0]), struct.pack( end+'c', ident[1]), struct.pack( end+'c', ident[2]), struct.pack( end+'c', ident[3]))
#       leads to an bug, because 4 chars do not equal an string of for chars:
#       fwrite.write( struct.pack( end+'s', ident)) #end+'4c'
       try:
           towrite = lambd(snap).astype('f').tostring()
           fwrite.write( struct.pack( end+'c', ident[0]) + struct.pack( end+'c', ident[1]) + struct.pack( end+'c', ident[2]) + struct.pack( end+'c', ident[3])) 
           fwrite.write( towrite  )#end+f
       except:
           print('save_snap:',ident,'not found')

       
   pairs  = [ ['HSML',lambda x:x.hsml],
#              ['MGAS',lambda x:x.mgas],
              ['U   ',lambda x:x.u],
              ['RHO ',lambda x:x.rho],
              ['ENDT',lambda x:x.endt],
              ['UUP ',lambda x:x.uup],
              ['UDOW',lambda x:x.udow],
              ['RUP ',lambda x:x.rup],
              ['RDOW',lambda x:x.rdow],
              ['DIVV',lambda x:x.divv],
              ['MACH',lambda x:x.mach],]
   for ident,lambd in pairs:
#       leads to an bug, because 4 chars do not equal an string of for chars:
#       fwrite.write( struct.pack( end+'s', ident)) #
       try:
           towrite = lambd(snap).astype('f').tostring()
           fwrite.write( struct.pack( end+'c', ident[0]) + struct.pack( end+'c', ident[1]) + struct.pack( end+'c', ident[2]) + struct.pack( end+'c', ident[3]))
           fwrite.write( towrite) 
       except:
           print('save_snap:',ident,'not found')

   fwrite.close()
   
   if verbose:
      print(' saving snapshot ... done')
     
   return 




class Snapshot:
   # A simple class for a snapshot of one resimulated Volume of MUSIC-2

   def __init__(self, name, radio_name='', transform=True, headerc='Mpc') :
       
      self.name       = name.replace('Additional2','') # Compability with old folder structure
      self.radio_name = radio_name
      self.head = None 
      self.hsml = None
      self.pos  = None 
      self.vel  = None
      self.mgas = None
      self.id   = None 
      self.u    = None
      self.rho  = None 
      self.endt = None 
      self.norm = None
      self.uup  = None
      self.udow = None 
      self.rup  = None
      self.rdow = None 
      self.divv = None 
      self.mach = None
      self.radi = None
      self.xray = None 
      self.transform  = transform
      self.headerc    = headerc

   def loaddata(self, term=False) : 
      if term: print('endian= ', '<') # args.endian
      load_snap ( self, ) 
      
   def savedata(self, savefolder, term=False):
      if term: print('endian= ', '<') # args.endian
      save_snap ( self, savefolder)    
      
# conversion factor gadget density to proper baryon density [cm-3]
def conversion_fact_gadget_rho_to_nb ( head ) :
   fact  = 1e10 * const.solar_mass / head['hubble'] / const.proton_mass
   fact /= ( 1e3 * const.parsec / head['hubble'] )**3.
   fact /= ( head['aexpan'] )**3. 
   return fact 
   
# conversion factor gadget density to proper baryon density [cm-3] as if the snapshot was taken at z=0
def conversion_fact_gadget_rho_to_nb_z0 ( head ) :
   fact  = 1e10 * const.solar_mass / head['hubble'] / const.proton_mass
   fact /= ( 1e3 * const.parsec / head['hubble'] )**3.
   #fact /= ( head['aexpan'] )**3. 
   return fact  
   
# electron density for fully ionized medium
# ......    Y = 4nHe / ( nH + 4nHe )
# ......    ne / nb = nH + 2nHe / ( nH + 4nHe )
# ......    nb / ntot = nH + 4nHe / ( 2nH + 3nHe )    (note ne = nH + 2nHe, fully ionized)
def conversion_fact_ne_per_nb( ) :	
	nHe_per_nH  = const.Y_helium / ( 4.0 * ( 1.0 - const.Y_helium ))
	ne_per_nb   = 1.0 + 2.0 * nHe_per_nH / ( 1 + 4.0 * nHe_per_nH )
	return ne_per_nb
   
def conversion_fact_nb_per_ntot ( ) :
	nHe_per_nH  = const.Y_helium / ( 4.0 * ( 1.0 - const.Y_helium ))
	nb_per_ntot = ( 1.0 + 4.0*nHe_per_nH ) / ( 2.0 + 3.0 * nHe_per_nH )
	return nb_per_ntot
   
# convert gadget specific energy density (U) to temperature [keV]
# ......   GADGET CONVERSION FROM U TO T
# ......   XH      = Zmgas(6,*)/mass_gas(*)
# ......   yHelium = (1. - XH)/(4.*XH)
# ......   mu      = (1 + 4.* yHelium)/ (1.+ yHelium + ne1gas )
# ......   temp    = GAMMA_MINUS1 * ugas * mu * 1.6726 / 1.3806 * 1.e-8 ; / BOLTZMANN  * PROTONMASS
# ......   temp    = temp * 1e10 ; UnitEnergy_in_cgs/UnitMass_in_g (to get T in Kelvin)
def conversion_fact_gadget_U_to_keV( head ) :
   nb_per_ntot = conversion_fact_nb_per_ntot( )
   fact  = const.kilometer_per_sec**2. * ( head['gamma'] - 1. ) * nb_per_ntot * const.proton_mass
   fact /= 1e3*const.electron_volt
   return fact


def comH_to_phys( head, z=None ) :
 
    if z is not None:
        return 1/(1+z)/head['hubble'] 
    else:
        return head['aexpan']/head['hubble']



