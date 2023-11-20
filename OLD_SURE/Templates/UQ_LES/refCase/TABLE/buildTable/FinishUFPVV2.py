#!/usr/bin/env python2.7
"""
This tools finishes the set up of a 4D UFPV table by building the ChiOverChist variable. (Reference needed)
It also does some miscellaneous changes in the names of certain parameters or fields.
It also adds a dummy flamelet at he second critical point of the SCurve, with the same values as the penultimate one, but the source terms are all naught
"""

from scipy.special import erfcinv, erfinv, gamma
from scipy import integrate
from math import exp, pi
import argparse
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib as mpl
font = {'family' : 'sans',
    'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)
import h5py


##################
#    Parsing     #
##################
parser = argparse.ArgumentParser(description="Finish UFPV tabulation for AVBP")
parser.add_argument('-t', '--table', required=True, help = "Path to table")
parser.add_argument('--plot'   , action='store_true', help = "plot Frond")
parser.add_argument('--noWrite', action='store_true', help = "Prevent effective encryption on table file")
parser.add_argument('--dummySR', required=True, help = "SR at which the flame cannot ignite. A dummy flamelet is created at this SR with same values except source terms which are set to 0")
args = parser.parse_args()

table     = args.table
plot      = args.plot
noWrite   = args.noWrite
dummySR   = float(args.dummySR)


##################
#   Parameters   #
##################
S_ZMin = 0.05


##################
#   Functions    #
##################
def Flam(Z, Zst):
  res = exp( 2.0 * (( erfcinv( 2.0*Zst ) )**2 - (erfcinv( 2.0*Z ) )**2 ) )
  return res
  
def FStraight(Z):
  res = 1.0 * exp( -2.0 * ( erfinv(2.0*Z - 1.0) )**2 )
  res/= pi
  return res

def integrand(zstar, a, b):
  return FStraight(zstar) * zstar**(a-1.0) * (1.0-zstar)**(b-1.0)

def FrondTh(Z, S_Z, Zst):
  a = Z * (1.0/S_Z - 1.0)
  b = a * (1.0/Z - 1.0)
  integral, error = integrate.quad(integrand,0,1, args=(a,b))
  res = gamma(a+b) * integral / ( FStraight(Zst) * gamma(a) * gamma(b) )
  return res

def linearInterp(x, x0, x1, y0, y1):
  return y0 + (y1-y0) * (x-x0) / (x1-x0)

def Frond(Z, S_Z, Zst):
  if Z <= 1e-6 or Z > (1.0 - 1e-6):
    return 0.0
  if S_Z > (1.0 - 1e-6):
    return 0.0
  if S_Z < 1e-6:
    return Flam(Z,Zst)
  if S_Z > S_ZMin:
    return FrondTh(Z, S_Z, Zst)
  else:
    y0 = Flam(Z,Zst)
    y1 = FrondTh(Z,S_ZMin, Zst)
    return linearInterp(S_Z, 0.0, S_ZMin, y0, y1)

def ComputeChist(SR, Zst):
  return SR / pi * exp(-2.0 * (erfcinv(2.0*Zst))**2)



##################
#   Computing    #
##################

f    = h5py.File(table, 'r+')
"""
Here, we assume that the variables of the table are:
1. Chist
2. Z
3. C
4. S_Z
"""
xChist  = f["Parameters/DIM1_INDEX"][:]
xZ      = f["Parameters/DIM2_INDEX"][:]
xC      = f["Parameters/DIM3_INDEX"][:]
xS_Z    = f["Parameters/DIM4_INDEX"][:]
zst     = float(f["Header/Z_st"][0])
nChist  = len(xChist)
nZ      = len(xZ)
nC      = len(xC)
nS_Z    = len(xS_Z)
ntot    = len(xChist) * len(xZ) * len(xC) * len(xS_Z)


# Building xChiOverChist
X, Y = np.meshgrid(xZ, xS_Z)
Res = np.zeros_like(X)
for iz in range(len(xZ)):
  for isz in range(len(xS_Z)):
    Res[isz,iz] = Frond(xZ[iz], xS_Z[isz], zst)

if not noWrite:
  # Miscellaneous changes:
  #f["Data/SIGMA_Z"]                        = f["Data/SIGMA_Z_tild"] 
  #f["Data/CHI_Z"]                          = f["Data/CHI_Z_tild"] 
  f["Data/DYK_DC"]                          = f["Data/DYK_DPROG_VAR"] 
  #f.create_dataset("Data/DYK_DZ", data=[0]*ntot)
  #del f["Data/SIGMA_Z_tild"] 
  #del f["Data/CHI_Z_tild"] 
  del f["Data/DYK_DPROG_VAR"] 
  f["Parameters/DIM3_NAME/"][...] = "C"

  # Writing ChiOverChist:
  myList = []
  for isz in range(len(xS_Z)):
    for c in xC:
      for iz in range(len(xZ)):
        for chist in xChist:
          myList.append(Res[isz,iz])
  f.create_dataset("Data/ChiOverChist", data=myList)
  ######################################################################
  # Adding a dummy flamelet at no-ignition SR with same properties as the last one except source terms which are set to 0
  ######################################################################
  # Changing chist indices :
  del f["Parameters/DIM1_INDEX"]
  newChist = ComputeChist(dummySR, zst)
  f.create_dataset("Parameters/DIM1_INDEX", data = np.append(xChist, newChist))
  del f["Parameters/NPTS1"]
  f.create_dataset("Parameters/NPTS1", data = [nChist+1])

  # Changing the dataset for every variable
  DataDict = {}
  def AddDummyFlameletInfo(name, obj):
      print(name)
      oriVect   = f["Data/"+ name][:]
      oriMat    = np.reshape(oriVect, (nChist, nZ, nC, nS_Z), order='F')
      newMat    = np.zeros((nChist+1, nZ, nC, nS_Z))
      newMat[:-1,:,:,:] = oriMat
      if 'W_' in name:
        pass
      elif name=='Chi_st':
        newMat[-1,:,:,:] = newChist
      else:
        newMat[-1,:,:,:] = oriMat[-1,:,:,:]
      newVect = np.reshape(newMat, newMat.size, order='F')

      DataDict[name] = newVect

  f["Data"].visititems(AddDummyFlameletInfo)

  for name in DataDict.keys():
    del f["Data/" + name]
    f.create_dataset("Data/" + name, data = DataDict[name])
 

f.close()

if plot:
  fig = plt.figure()
  ax = plt.axes(projection='3d')
  ax.plot_wireframe(X, Y, Res, color='black')
  ax.set_title(r'$\frac{\chi}{\chi_{st}}$')
  ax.set_xlabel(r'$\tilde{Z}$', labelpad=20)
  ax.set_ylabel(r'$S_Z$', labelpad=20)
  ax.set_zlabel(r'$\frac{\chi}{\chi_{st}}$', labelpad=20)
  plt.show()











