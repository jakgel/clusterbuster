import numpy as np
import matplotlib.mlab as mlab
import math
import random


def sigmoid(x):
  return 1 / (1 + np.exp(-x))


def interpolate(x,y,z,subgrid=8):
  
  # Normalize coordinate system
  def normalize_x(data):
      data = data.astype(np.float)
      return (data - xmin) / (xmax - xmin)

  def normalize_y(data):
      data = data.astype(np.float)
      return (data - ymin) / (ymax - ymin)
  
  xmin, xmax = x[0], x[-1]#0, 1
  ymin, ymax = y[0], y[-1]#0, 10
  
  # Size of regular grid
  ny, nx = len(y)*subgrid, len(x)*subgrid
  xnew   = np.repeat(x, len(y)) # Kind of 1 1 1 2 2 2 3 3 3 #np.tile(x,   len(y))    # Kind of 1 2 3 1 2 3 1 2 3
  ynew   = np.tile(y,   len(x)) # Kind of 1 2 3 1 2 3 1 2 3 #np.repeat(y, len(x))    # Kind of 1 1 1 2 2 2 3 3 3
  znew   = z.flatten()
  

  # Generate a regular grid to interpolate the data.
  xi = np.linspace(xmin, xmax, nx)
  yi = np.linspace(ymin, ymax, ny)
  xii, yii = np.meshgrid(xi, yi)
  
  # Interpolate using delaunay triangularization 
  zi = mlab.griddata(xnew,ynew,znew,xii,yii,interp='linear') # without 'linear' it would be like that
    
  return xi, yi, zi
  
  
  
# from http://stackoverflow.com/questions/8914491/finding-the-nearest-value-and-return-the-index-of-array-in-python
def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx
    
    

# Definition from english wikipedia for rotation matrixes following the right hand rule - which codifies their alternating signs
def RotM_x(a): 
   # 3D Rotation Matrix, of what kind?
   ca = np.cos(a)
   sa = np.sin(a)
   
   return np.array( [ [ 1,   0,  0], 
                      [ 0,  ca, sa], 
                      [ 0, -sa, ca] ] )
   
def RotM_y(a): 
   # 3D Rotation Matrix, of what kind?
   ca = np.cos(a)
   sa = np.sin(a)
   
   return np.array( [ [ ca, 0, -sa], 
                      [  0, 1,   0], 
                      [ sa, 0,  ca]] )

def RotM_z(a): 
   # 3D Rotation Matrix, of what kind?
   ca = np.cos(a)
   sa = np.sin(a)
   
   return np.array( [ [ ca, sa, 0], 
                      [-sa, ca, 0], 
                      [ 0,   0, 1] ] )

def Kep3D(Omega,i,omega): 
   # Kepleria Rotation Matrix, To rotate an vector from reference direction into the eulerian coordinate system
   
   # My defintion has the negative sin values, compared to sebstian onces, which is why I should use exactly the negitive angular values
   return  (RotM_x(Omega).dot(RotM_y(i))).dot(RotM_z(omega))

def Kep3D_seb(Omega,i,omega):  #Follows the 'normal' convention
  return Kep3D(-Omega,-i,-omega)

def findAnglesBetweenTwoVectors3D(v1s, v2s):   # taken from a python forum http://codereview.stackexchange.com/questions/54347/calculating-element-wise-the-angles-between-two-lists-of-vectors
    dot = np.einsum('ijk,ijk->ij',[v1s,v1s,v2s],[v2s,v1s,v2s])
    return np.arccos(dot[0,:]/(np.sqrt(dot[1,:])*np.sqrt(dot[2,:])))
    
    
def findAnglesBetweenTwoVectors2D(v1, v2):   # taken from a python forum
    #a  = np.arctan2(v1s[0] ,v1s[1])
    #aa = np.arctan2(-v1s[0],-v1s[1])
    #b  = np.arctan2(v2s[0],v2s[1])
    
    #mhhh  = [np.absolute(a-b),np.absolute(aa-b)]
    #return min(mhhh)
    a_v1 = np.arctan2( v1[1],v1[0] )
    a_v2 = np.arctan2( v2[1],v2[0] )
    diff  = a_v1-a_v2
    return      np.max( [diff,2*np.pi-diff] )   #np.min( np.absolute( [ (a_v1-a_v2), (a_v1+a_v2)] ) )
 
 
 
#taken from http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
def returnangle(u, v):
#Returns the angle beetween two vectors  in [0,180] degree
    c     = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v) # -> cosine of the angle
    angle = np.arccos(np.clip(c, -1, 1)) 
    return np.degrees(angle) 
    
 
def returnangle_quarter(u, v): 
    return       min( returnangle(u,v),returnangle(u,-v) )
    
    
def angle2vec(angle, deg=True): 
# This folows the astronomical coordinate system convention with an clockwise polar coordinate system starting with theta=0 for [0,1] an vector to the north
  
  if deg:
    ang = np.radians(angle) 
  else:
    ang = angle

  return [np.sin(ang),np.cos(ang)] #np.asarray([np.cos(ang),np.sin(ang)]) #np.column_stack((np.cos(ang),np.sin(ang))) 


def MinAngle( th1, th2):
    diff  = th1-th2
    return   np.max( np.concatenate( (diff,360-diff), axis = 0 )  ) #np.min( np.absolute( [ (a_v1-a_v2), (a_v1+a_v2)] ) )
    
def MinAngle_quarter( th1, th2):
    diff  = th1-th2
    return   np.min(np.concatenate( ( [(180+diff) % 180], [(180-diff) % 180] ), axis=0 ), axis=0 )
  
    
def polar_from_to(array,borders):
    
        start,end = borders
        N = array.shape[0]
        

#        array = np.asarray( range(array.shape[0]-1) )  # Debugging
        
        if start >= 0 and end < N:
            array_polar = array[borders[0]:borders[1]]
        elif start < 0 and end < N:
            array_polar = np.append(  array[N+start-1:N] , array[0:end+1] )     
        elif start >= 0 and end >= N:
            array_polar = np.append(  array[start:N] , array[0:(end % N)+1] )   
        else:
            print('Something went wrong in polar_from_to(array,%i,%i)' % (start,end))
        
#        print(start, end, array_polar) # Debugging
        return array_polar


def compute_moments(sparseA, sparseD, sparseW, xy_swapped = True):

    if len(sparseD) > 0:
        sparseD = [1 for a in sparseA] # DEBUGGING

        if xy_swapped:
            x = np.cos(sparseA)*sparseD # * -1 next test
            y = np.sin(sparseA)*sparseD
        else:
            x = np.sin(sparseA)*sparseD
            y = np.cos(sparseA)*sparseD

        moment_00 = np.sum(sparseW)
        moment_10 = np.sum(sparseW*x)  # x*x
        moment_01 = np.sum(sparseW*y)  # y*y
        moment_11 = np.sum(sparseW*x*y)
        moment_20 = np.sum(sparseW*x*x)
        moment_02 = np.sum(sparseW*y*y)

        _xm = moment_10/moment_00
        _ym = moment_01/moment_00

        moment_central_11 = moment_11/moment_00 - _xm*_ym
        moment_central_20 = moment_20/moment_00 - _xm*_xm
        moment_central_02 = moment_02/moment_00 - _ym*_ym
        moment_angle = .5*np.arctan2(2*moment_central_11, moment_central_20-moment_central_02)
        #moment_angle = np.arctan(moment_central_11,moment_central_20-moment_central_02)
        #print('Moments:', moment_00, moment_10, moment_01, moment_11, moment_20, moment_02)
        #print('Central:', moment_central_11, moment_central_20, moment_central_02)
        #print('Derived', moment_angle)

        #if xy_swapped:
        #    moment_angle + np.pi/2

    else:
        moment_angle = 0


    return moment_angle




def nextTime(rateParameter):
    return -math.log(1.0 - random.random()) / rateParameter



#=== TODO: Some cosmological stuff which could be put somewhere else ===#
def Mvir2Lx(MassBins,h=0.7): #Made for logarithmic values
   # Using the scalling relation of Reichert et al for z=0.0
  Lx    = 1.52*MassBins+1.52*np.log10(0.54*1e-15/h)+np.log10(0.88*2.5*1e45)    #W/Hz read file ... convert to masses ... convert to erg/s 
  return Lx
  
  



