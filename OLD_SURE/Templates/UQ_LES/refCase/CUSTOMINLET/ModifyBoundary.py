#!/usr/bin/env python2.7

import numpy as np
import h5py
import os
import sys
sureHome = os.environ['SURE_HOME']
sys.path.insert(0, sureHome)
from SurrogateModelling.OneDSpline import OneDSpline

#import matplotlib.pyplot                                  as plt
#import matplotlib                                         as mpl
#font = {'family' : 'sans',
#    'serif': 'helvet',
#    'weight' : 'bold',
#    'size'   : 20}
#mpl.rc('font', **font)
#mpl.rc('text', usetex=True)


def importFile(fileName):
  delimiterIn = ' '
  nHeaderLines = 4
  X = []
  Y = []
  with open(fileName, 'r') as f:
    content = f.readlines()
    for lineIdx in range(len(content)-1):
      if lineIdx >= nHeaderLines:
        line = content[lineIdx]
        row = line.split()
        X.append(float(row[0]))
        Y.append(float(row[1]))
  return np.array(X), np.array(Y)


# ----- Define work parameters -----
meshFile    = "../MESH/mesh.mesh.h5"
boundFile   = "../SOLUTBOUND/mesh.solutBound.h5"
UFile       = 'velocityProfile'
intensFile  = 'intensityProfile'
patchName   = 'inlet_fuel'
patchNumber = 1

# ----- Read velocity and intensity profiles -----
z1, U = importFile(UFile)
z2, I = importFile(intensFile)
if len(z1) == len(z2):
  Urms = U * I / 100

# ----- Build surrogates for U and Urms -----
USpline     = OneDSpline(ID='USpline', verbosity=0)
UrmsSpline  = OneDSpline(ID='UrmsSpline', verbosity=0)
USpline.Build(samples=z1, evaluations=U)
UrmsSpline.Build(samples=z1, evaluations=Urms)


newZ    = np.linspace(z1.min(), z1.max(), 1000)
newU    = USpline.Evaluate(newZ)
newUrms = UrmsSpline.Evaluate(newZ)

# ----- Plot profiles -----
#fig = plt.figure(1, figsize=(14,8))
#ax1 = fig.add_subplot(121)
#ax1.scatter(z1, U, color='gray', marker='x', label='Original') 
#ax1.plot(newZ,newU, color='k', label='surrogate')
#ax1.set_xlabel(r'$r$')
#ax1.set_ylabel(r'$U$')
#ax1.set_xlim([z1.min(), z1.max()])
#ax1.legend(loc='best')
#plt.legend(loc='best')
#for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
#  item.set_fontsize(18)
#ax2 = fig.add_subplot(122)
#ax2.scatter(z2, U*I/100, color='gray', marker='x', label='Original') 
#ax2.plot(newZ,newUrms, color='k', label='surrogate')
#ax2.set_xlabel(r'$r$')
#ax2.set_ylabel(r'$U_{rms}$')
#ax2.set_xlim([z2.min(), z2.max()])
#ax2.legend(loc='best')
#plt.legend(loc='best')
#for item in [ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels():
#  item.set_fontsize(18)
#plt.show()
#plt.cla()
#plt.clf()
#plt.close()

# ----- Read coordinates -----
mesh		= h5py.File(meshFile, 'r')
xx		  = mesh["Patch"][str(patchNumber)]['Coordinates/x'][:]
yy		  = mesh["Patch"][str(patchNumber)]['Coordinates/y'][:]
zz		  = mesh["Patch"][str(patchNumber)]['Coordinates/z'][:]
mesh.close()  

# ----- Compute values -----
rList = np.sqrt(yy*yy+zz*zz)
uList = USpline.Evaluate(rList)
urmsList = UrmsSpline.Evaluate(rList)

# ----- Write solutbound -----
boundPath = 'Patch_00' + str(patchNumber) + '-' + patchName
bound	= h5py.File(boundFile, 'r+')
del bound[boundPath]['U']
del bound[boundPath]['Urms']
bound.create_dataset(boundPath + '/U', data=uList)
bound.create_dataset(boundPath + '/Urms', data=urmsList)
bound.close()  
