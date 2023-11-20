#!/usr/bin/env python2.7


import numpy as np
import matplotlib.pyplot as plt
import Link
import openturns as ot
from Sampling.MC   import *
from Numerics.Algebra import Orthonormalize
from Numerics.Algebra import CompleteBasis
from Numerics.Algebra import ChangeBasis
from Numerics.Algebra import ComputeActiveVectors
from Numerics.Algebra import VectorProjection
from Numerics.Algebra import norm2
from DimensionReduction.ActiveDirection import *
from DimensionReduction.Preprocessing import *
from SurrogateModelling.PCE import *
from PoolUsefulScripts.OK_model import *
import matplotlib as mpl
font = {'family' : 'sans',
   'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)



#------------------------------------------------------------ Meta-parameters
#nSample      = 50
nSample      = 100
firstSample  = 0
nDim         = 34
distribution = ot.ComposedDistribution([ot.Uniform(0.0, 1.0)] * nDim)
myThreshold  = 1.0e-2
#myThreshold  = 1.0e-1
#myThreshold  = -1.0

#------------------------------------------------------------ Miscellaneous plot functions
def ReductionPlot1D(x,y, sur=None):
  fig = plt.figure(1, figsize=(7,5))
  ax1 = fig.add_subplot(111)
  ax1.scatter(x, y, color='k', alpha = 0.5, marker=',', s=1, label='original data')
  ax1.set_xlabel("reduced variable")
  ax1.set_ylabel("data")
  ax1.set_title("Reduction plot")
  ax1.legend(loc='best')

  if sur is not None:
    nReSample = 1000
    x1 = np.linspace(x.min(), x.max(), nReSample)
    y1 = np.zeros(nReSample)
    for idx in range(nReSample):
      y1[idx] = sur.predict(x1[idx])
    ax1.plot(x1,y1, label='Surrogate')

  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)
  plt.show()
  plt.close()
#
def ReductionPlot2D(x,y,z, sur=None):
  #from matplotlib import animation
  from mpl_toolkits.mplot3d import Axes3D
  fig = plt.figure()
  ax = Axes3D(fig)
  ax.scatter(x, y, z, marker='3', s=20, c="k", alpha=0.6, label='Original data projected on the 2 directions')

  if sur is not None:
    nReSample = 100
    #x1 = np.linspace(x.min(), x.max(), nReSample)
    #y1 = np.linspace(y.min(), y.max(), nReSample)
    x1 = np.linspace(-2, 2, nReSample)
    y1 = np.linspace(-2, 2, nReSample)
    z1 = np.zeros((nReSample, nReSample))
    X,Y = np.meshgrid(x1,y1)
    for idx1 in range(nReSample):
      for idx2 in range(nReSample):
        z1[idx1,idx2] = sur.predict([x1[idx1],y1[idx2]])
    ax.plot_wireframe(Y,X, z1, label='Surrogate')

  plt.legend(loc='best')
  plt.show()
  plt.close()


#------------------------------------------------------------ Load samples
resultFile = 'results/samples_fusion.h5'
result = h5py.File(resultFile, 'r')
o0 = np.log(result["Data/QoI"][:nSample])
s0 = result["Data/Inputs"][:nSample, :].T
result.close()
PostProc = CenterScaling('myPostProc', o0, 'AUTO')
obs = PostProc.ComputePro(o0)
PreProc = CenterScaling('myPreProc', s0, 'AUTO')
sam = PreProc.ComputePro(s0).T

####################################################################################################
##   Full MAD baby !   ###
##########################

#------------------------------------------------------------ Initialization
counter = 0
xCurrentInput     = None
xRemInput         = np.zeros((nSample, nDim))
xRemInput[:,:]    = sam[:,:]
xAD               = []
xRemObs           = np.zeros_like(obs)
xRemObs[:]        = obs[:]
xBasis            = []
xxActiveVariable  = []
xKrig             = []


#------------------------------------------------------------ Loop
while counter < nDim:
  transferMatrix = []
  remainingDim = nDim - counter

  #---------------------------------------------------------- Compute current active direction
  #prt(' Finding active direction ' + str(counter) +'...', 'None', True)
  ID = 'AD_' + str(counter)
  AD = ActiveDirection(ID=ID, totalDim=remainingDim, verbosity=-1)
  AD.Build(xRemInput, xRemObs)
  xAD.append(AD)

  #---------------------------------------------------------- Orthonormalizing the basis
  #prt(' Orthonormalizing the basis...', 'None', True)
  transferMatrix.append(xAD[-1].xWeight)
  completeBasis = Orthonormalize(CompleteBasis(np.array(transferMatrix)))
  # No need to update the active direction as neither the completion nor the orthonormalization change the first vector of the basis
  xBasis.append(completeBasis)

  #---------------------------------------------------------- Computing the coordinates of the remaining inputs in the new basis
  xRemInputNB = np.zeros_like(xRemInput)
  for i in range(nSample):
    xRemInputNB[i,:] = ChangeBasis(xRemInput[i,:], completeBasis)

  #---------------------------------------------------------- Extracting the active variables and the remaining input at this stage
  xAV = xRemInputNB[:,0]
  xxActiveVariable.append(xAV)
  xRemInput = None
  xRemInput = np.zeros((nSample, nDim-(counter+1)))
  xRemInput[:,:] = xRemInputNB[:,1:]

  #---------------------------------------------------------- Computing the active vectors at this stage
  xProjectedInput = ComputeActiveVectors(xxActiveVariable, xBasis)

  #---------------------------------------------------------- Building surrogate
  #currentKrig = OK_model("Matern52", np.array(xxActiveVariable).T, obs, givenTheta=[1e-12]*(counter+1), nugget=1.0) # Works nicely
  #currentKrig = OK_model("Matern52", np.array(xxActiveVariable).T, obs, givenTheta=[1e-12]*(counter+1), nugget=1e-10)
  #currentKrig = OK_model("Matern52", np.array(xxActiveVariable).T, obs, givenTheta=[1e-50]*(counter+1))
  #currentKrig = OK_model("Matern52", np.array(xxActiveVariable).T, obs)
  currentKrig = OK_model("Matern52", np.array(xxActiveVariable).T, obs, optiStrat="GeomDichotomy", optiTarget="NRMSE", nugget = 1.0)
  #currentKrig = OK_model("Matern52", np.array(xxActiveVariable).T, obs, nugget=0.01)
  print(currentKrig.theta)
  print('NRMSE = ' + str(currentKrig.computeNRMSE(currentKrig.theta)))
  xKrig.append(currentKrig)
  approx = np.zeros(nSample)
  for i in range(nSample):
    approx[i] = currentKrig.predict(np.array(xxActiveVariable).T[i,:])
  xRemObs[:] = obs[:] - approx[:]



  # ---------------------------------------------------------- Plots for 1D and 2D
  counter += 1
  if counter == 1:
    OneDabs = xxActiveVariable[0]
    ReductionPlot1D(OneDabs,obs, sur=xKrig[0])

  elif counter ==2:
    TwoDabs = [xxActiveVariable[0], xxActiveVariable[1]]
    ReductionPlot2D(TwoDabs[0], TwoDabs[1], obs, sur=xKrig[1])

  # ---------------------------------------------------------- Test for exiting the loop
  testRatio = np.var(xRemObs) / np.var(obs)
  print("Proportion of the variance discarded: %.3G" % (testRatio * 100) + '%')
  if testRatio < myThreshold:
    break


