#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:34:32 2017

@author: jakobg
"""

from __future__ import division,print_function

import numpy as np

class entry:
  
    def __init__(self, val, dic, label=False,):
        self.value = val       # value e.g. zero moment
        self.dic = str(dic)    # Name of this variable

        if label:
            self.label = str(label)
        else:
            self.label = self.dic      #label of this variable, dictionary value if not stated otherwise

    def __call__(self):
        return self.value

    
class norm:
    """ A class for a normalizing constant """
    def __init__(self, ntype='R200', Nexp=1):
        self.ntype = ntype
        self.Nexp = Nexp
    
    
class reference:
    
    """ This class was introduced to manage references """
    
    def __init__(self, ident, rtype='text', page=None, nr=None):

        self.ident  = ident   # Identifier in the citations.bib
        self.rtype  = rtype   # table or section
        self.nr     = nr      # section/table number

class measurand(entry):
    
    """ 
    An class that inherits from 'entry' and describes measures wit units and errors.
    Inspired by http://www.python-course.eu/python3_magic_methods.php it also provides some magic operation functionalities that should make it more usable like float numbers.
    """
    
    def __init__(self, value, dic, label=False, std=0, skew=0, vrange=[None,None], vmin=None, vmax=None, un='arbitrary unit', distr='gauss', ref=None):
        
        entry.__init__(self, value, dic, label)

        self.unit = un # A unit, here as a string

        #== variable range
        self.vrange = vrange
        self.set_vrange(vmin, vmax)

        #=== Errors
        self.distr = distr   # distribution type
        #self.distr = ['norm','linear','weird','lgnorm']
        
        self.set_std(std)   # first moment, can be one or two dimensional
        self.skew = skew    # second momnt
        self.ref = ref      #reference
    
    def set_vrange(self, vmin, vmax):
        if vmin is not None: self.vrange[0] = vmin
        if vmax is not None: self.vrange[1] = vmax               
    
    def set_std(self, std):
        try:
            self.std = [std[0],std[1]]
        except:
            self.std = [std, std]
                 
    def std_add(self, std):  # Also use reversed for other operations
        self.std[0] = np.sqrt(self.std[0]**2 + std[0])
        self.std[1] = np.sqrt(self.std[1]**2 + std[1])


    def unit_string(self):

        if self.unit is None:
            return ''
        else:
            return '[%s]' % (self.unit)

    def labels(self, log=False):
        """ Returns the labels of measurand objects  
            Use 'log=True' if you want to indicate that you use the logarithmic values
        """
        if self.unit is None:
            unit = ''
        else:
            unit = ' %s' % (self.unit_string())

        if log:
            return '$log_{10}($%s%s$)$' % (self.label, unit)
        else:
            return '%s%s' % (self.label, unit)
        
          
    def inrange(self, limits):
      return (self.value>limits[0]) & (self.value < limits[1])     
        
    def __str__(self):
        # return("The cluster %12s at dec=%6.2f and rec = %6.2f has an z of %5.3f" % (self.name, self.dec, self.rec, self.z))
        return '%.3e' % (self.value)


    def __add__(self, other):
        """
        Defines the '+' operator.
        If other is a measurand object the currency values 
        are added and the result will be the unit of 
        self. 
        """      
        if isinstance(other, measurand):
            return self.value + other.value
        else:
            return self.value + other 

    def __radd__(self, other):
        return self.__add__(other)    #measurand.__add__(self,other)

    def __iadd__(self, other):
        """
        Similar to __add__
        """
        if  isinstance(other, measurand):   
            self.value += other.value
        else:
            self.value += other
        return self
    
    def __sub__(self, other):
        """
        Defines the '-' operator.
        If other is a measurand object the currency values 
        are added and the result will be the unit of 
        self. 
        """
        if  isinstance(other, measurand):   
            return self.value - other.value
        else:
            return self.value - other
        
    def __rsub__(self, other):
        return self.__sub__(other)   #self.__sub__(self,other) 
    
    def __isub_(self, other):
        """
        Similar to __sub__
        """
        
        if  isinstance(other, measurand):   
            self.value -= other.value
        else:
            self.value -= other
        return self
    
    
    def __mul__(self, other):
        """
        Multiplication is only defined as a scalar multiplication.
        """
        if  isinstance(other, measurand):   
            return self.value * other.value
        else:
            return self.value * other 
        
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __imul__(self, other):
        if  isinstance(other, measurand): 
            self.value *= other .value
        else:
            self.value *= other
        return self
    

    def __truediv__(self, other):
        """
        Division is only defined as a scalar division.
        """
        if  isinstance(other, measurand): 
            return self.value / other.value
        else:
            return self.value / other
        
    def __rtruediv__(self, other):
        return self.__div__(other)
    
    def __itruediv__(self, other):
        if  isinstance(other, measurand): 
            self.value /= other.value
        else:
            self.value /= other
        return self
    

    def __div__(self, other):
        """
        Division is only defined as a scalar division.
        """
        if  isinstance(other, measurand): 
            return self.value / other.value
        else:
            return self.value / other
        
    def __rdiv__(self, other):
        return self.__div__(other)
    
    def __idiv__(self, other):
        if  isinstance(other, measurand): 
            self.value /= other.value
        else:
            self.value /= other
        return self
    
    def __pow__(self, other):
        """
        Power is only defined as a scalar division.
        """
        if  isinstance(other, measurand): 
            raise TypeError("unsupported operand type(s) for **: 'measurand' and " + type(other).__name__)   
        else:
            return self.value ** other

    def __neg__(self):
        return -self.value
    
    def __pos__(self):
        return self.value
    
    def __lt__(self, other): 
        if  isinstance(other, measurand):   
            return self.value > other.value
        else:
            return self.value < other

    def __gt__(self, other): 
        if  isinstance(other, measurand):   
            return self.value > other.value
        else:
            return self.value > other
        
    
    def __abs__(self):
        raise TypeError("unsupported operand type(s) for abs(): 'measurand'")  
        
    def __invert__(self):
        raise TypeError("unsupported operand type(s) for invert(): 'measurand'")  


    
class Histogram:
    
  def __init__(self, nbins=30, fromto = [0.,1.5]):


    self.nbins = nbins
    self.fromto = fromto
    self.width = np.abs(self.fromto[0]-self.fromto[1])/self.nbins
    self.bins  = np.linspace(self.fromto[0]             , self.fromto[1]           , self.nbins+1) 
    self.ticks = np.linspace(self.fromto[0]+self.width/2, self.fromto[1]-self.width, self.nbins+0) 
    
    self.hist = np.zeros( (nbins,) )
    
  def __str__(self):
    # return("The cluster %12s at dec=%6.2f and rec = %6.2f has an z of %5.3f" % (self.name, self.dec, self.rec, self.z))
    return 'Nothing specified in def__str__(self) for <Histogram>'
    #return("The cluster %12s at dec=%11s and ra = %9s has an z of %5.3f" % (self.name, self.dec, self.ra, self.z))      
    
  def _like(self):
      return Histogram(nbins=self.nbins, fromto= self.fromto )
    
         
class HistogramDD:
  """ By default it is an 3D Histogram """
  # rather in die format dimensions x , each dimensions has a from, to and a number of bins, as same as the widh  
  def __init__(self, Ndims=None, nbins=[30,30,30], fromto = [[0.,1.5]]*3, axisM=( measurand(0.,''),measurand(0.,''),measurand(0.,'') ), norm=norm() ):
     
    if Ndims is not None:
        nbins = [30 for i in np.arange(Ndims)]
    else:
        Ndims = len(nbins)
          
    self.nbins  = nbins
    self.fromto = fromto
    self.width  = [np.abs(self.fromto[i][0]-fromto[i][1])/nbins[i] for i in np.arange(Ndims)]
    self.bins   = [np.linspace(self.fromto[i][0]                , self.fromto[i][1]              , nbins[i]+1) for i in np.arange(Ndims)]
    self.ticks  = [np.linspace(self.fromto[i][0]+self.width[i]/2, self.fromto[i][1]-self.width[i], nbins[i]+0) for i in np.arange(Ndims)] 
    self.axisM  = axisM   # An 3 dimensional array of measurands
     
    self.hist = np.zeros(nbins)
    self.norm = norm
    
  def __str__(self):
    # return("The cluster %12s at dec=%6.2f and rec = %6.2f has an z of %5.3f" % (self.name, self.dec, self.rec, self.z))
    return 'Nothing specified in def__str__(self) for <Histogramm>'
    #return("The cluster %12s at dec=%11s and ra = %9s has an z of %5.3f" % (self.name, self.dec, self.ra, self.z))     
    
  def _like(self):
      return HistogramDD(nbins=self.nbins, fromto= self.fromto, norm=self.norm )
    
class Histogram2D(HistogramDD):
      
  # rather in die format dimensions x , each dimensions has a from, to and a number of bins, as same as the widh  
  def __init__(self, nbins=[30,30], fromto = [[0.,1.5]]*2, axisM=( measurand(0.,''),measurand(0.,'')), norm=norm() ) :
    HistogramDD.__init__(self, nbins=nbins, fromto=fromto, axisM=axisM, norm=norm)   
  
  
class Histogram3D(HistogramDD):    
    
  # rather in die format dimensions x , each dimensions has a from, to and a number of bins, as same as the widh  
  def __init__(self, nbins=[30,30,30], fromto = [[0.,1.5]]*3, axisM=( measurand(0.,''),measurand(0.,''),measurand(0.,'') ), norm=norm() ):
     
    HistogramDD.__init__(self, nbins=nbins, fromto=fromto, axisM=axisM,norm=norm)   