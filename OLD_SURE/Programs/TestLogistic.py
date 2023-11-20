#!/usr/bin/env python2.7

from __future__ import print_function
import h5py
import time
import Link
import numpy                                      as np
import openturns                                  as ot
import matplotlib.pyplot                          as plt
from SurrogateModelling.Logistic      import *
import Utils.PrettyPlots                          as pp
import matplotlib as mpl
font = {'family' : 'sans',
    'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
t0 = time.time()



counter = 0

#############################################################
# STEP 0: Retrieving the sampling data
#############################################################
resultFile = 'testData/solutFile.h5'
result  = h5py.File(resultFile, 'r')
xEval = np.log(result["Data/QoI"][:])
xSample = result["Data/Inputs"][:,:]
result.close()
nSample = len(xEval)

#############################################################
# STEP 1: Building a binary dataset
#############################################################
xBinaryEval = np.zeros_like(xEval)
mean = np.mean(xEval)
for i in range(nSample):
  if xEval[i] > mean:
    xBinaryEval[i] = 1
  else:
    xBinaryEval[i] = 0


#############################################################
# STEP 2: Building a logistic surrogate
#############################################################
print("Initialising Logistic")
myDummyLogistic = Logistic(ID="toto", verbosity=1)
print("Building Logistic")
myDummyLogistic.Build(xSample, xBinaryEval)
print("Logistic built")


##############################################################
# STEP 3: Evaluating the surrogate
#############################################################
print("Evaluating Logistic")
xReEval       = myDummyLogistic.Evaluate(xSample)
xReProbaEval  = myDummyLogistic.Evaluate(xSample, proba=True)
print("Logistic evaluated")

##############################################################
# STEP 4: Some plots
#############################################################
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(12,8))
ax = Axes3D(fig)
ax.scatter(xSample[:,0], xSample[:,1], xBinaryEval[:],      marker=',', s=1, c="k", alpha=0.7, label='Original')
ax.scatter(xSample[:,0], xSample[:,1], xReEval[:],          marker=',', s=1, c="r", alpha=0.7, label='Predicted')
ax.scatter(xSample[:,0], xSample[:,1], xReProbaEval[:, 1],  marker=',', s=1, c="g", alpha=0.9, label='Predicted probability')
#ax.plot_trisurf(xSample[:,15], xSample[:,16], xReProbaEval[:, 1])
ax.set_xlabel(r'$\xi_{16}$', labelpad=15)
ax.set_ylabel(r'$\hat{T}$',  labelpad=15)
ax.set_zlabel(r'$\delta$',   labelpad=15)
plt.title('Logistic surrogate test')
plt.legend(loc='best')
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.view_init(elev=7, azim=-47)
plt.show()
