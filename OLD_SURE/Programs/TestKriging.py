#!/usr/bin/env python2.7

from __future__ import print_function
import h5py
import time
import Link
import numpy                                      as np
import openturns                                  as ot
import matplotlib.pyplot                          as plt
from SurrogateModelling.Kriging      import *
from SurrogateModelling.OneDSpline      import *
from DimensionReduction.ActiveDirection import *
import Utils.PrettyPlots                          as pp

t0 = time.time()

#nSample = 5000
nSample = 2000


counter = 0

for T in [913.5]:
  #############################################################
  # STEP -1: Retrieving the sampling data
  #############################################################
  tFolderName	= './testData/T_' + str(T)
  print(tFolderName)
  resultFile = tFolderName + '/result.h5'
  result  = h5py.File(resultFile, 'r')
  xEval = np.log(result["Data/QoI"][:nSample])
  xSample = result["Data/NormalInputVars"][:nSample,:]
  result.close()
  
  ##Uncertain dimension and size of the dataset
  #dim = len(points[0,:])
  #dataSize = len(evaluations)
  
  
  #############################################################
  # STEP 0: Creating the Active Direction object
  #############################################################
  objAD = ActiveDirection(ID='tata', totalDim=xSample.shape[1], verbosity=1)
  objAD.Build(xSample, xEval)
  objAD.Save()

  #############################################################
  # STEP 1: Constructing the list of reduced samples
  #############################################################
  xRedSample = objAD.ComputeReducedFromTotal(xSample)

  #############################################################
  # STEP 2: Building a spline surrogate from the reduced var
  #############################################################
  objSpl = OneDSpline(ID='toto', verbosity=1)
  objSpl.Build(samples=xRedSample, evaluations=xEval)

  #############################################################
  # STEP 3: Building a Kriging surrogate from the reduced var
  #############################################################
  evals = xEval
  samples = np.zeros((nSample,1))
  for i in range(nSample):
    samples[i, 0] = xRedSample[i]

  objKrig = Kriging(ID='titi', verbosity=1, distribution=ot.ComposedDistribution([ot.Normal(0.0,1.0)]))
  objKrig.verbosity = True
  #objKrig.DefineTrend('Constant')
  objKrig.DefineTrend('Linear')
  #objKrig.DefineTrend('Quadratic')
  objKrig.DefineCovariance('Matern')
  objKrig.SetCovarianceOptimization()
  #objKrig.SetNoiseOptimization()
  objKrig.DefineNoise(1e-1)
  objKrig.Build(samples, evals)
  xKrigEval = np.zeros(nSample)
  for i in range(nSample):
    xKrigEval[i] = objKrig.Evaluate([xRedSample[i]])[0]


  #############################################################
  # STEP 4: Some plots
  #############################################################
  x1 = np.linspace(xEval.min(), xEval.max(), 1000)

  fig = plt.figure(1, figsize=(11,7))
  pp.AdjustSubplotSpacings()

  ax1 = fig.add_subplot(221)
  ax1.scatter(xEval, objSpl.Evaluate(xRedSample), color='k', alpha = 1.0, marker=',', s=1, label='Spline surrogate') 
  ax1.plot(x1,x1, color='g', label='Ideal')
  ax1.set_xlabel("Model evaluation")
  ax1.set_ylabel("Surrogate evaluation")
  ax1.set_title("Summary plot (Spline)")
  ax1.legend(loc='best')
  pp.CorrectAxSize(ax1)

  ax2 = fig.add_subplot(222)
  ax2.scatter(xEval, xKrigEval, color='k', alpha = 1.0, marker=',', s=1, label='Kriging Surrogate')
  ax2.plot(x1,x1, color='g', label='Ideal')
  ax2.set_xlabel("Model evaluation")
  ax2.set_ylabel("Surrogate evaluation")
  ax2.set_title("Summary plot (Kriging)")
  ax2.legend(loc='best')
  pp.CorrectAxSize(ax2)


  x2 = np.linspace(xRedSample.min(), xRedSample.max(), 1000)
  xKrigEval2 = np.zeros_like(x2)
  for i in range(len(x2)):
    xKrigEval2[i] = objKrig.Evaluate([x2[i]])[0]


  ax3 = fig.add_subplot(223)
  ax3.scatter(xRedSample, xEval, color='k', alpha = 1.0, marker=',', s=1, label='Samples') 
  ax3.plot(x2, objSpl.Evaluate(x2), color='g', label='Spline surrogate')
  ax3.set_xlabel("Reduced variable")
  ax3.set_ylabel("QoI")
  ax3.set_title("Model plot (Spline)")
  ax3.legend(loc='best')
  pp.CorrectAxSize(ax3)


  ax4 = fig.add_subplot(224)
  ax4.scatter(xRedSample, xEval, color='k', alpha = 1.0, marker=',', s=1, label='Samples') 
  ax4.plot(x2, xKrigEval2, color='g', label='Kringing surrogate')
  ax4.set_xlabel("Reduced variable")
  ax4.set_ylabel("QoI")
  ax4.set_title("Model plot (Kriging)")
  ax4.legend(loc='best')
  pp.CorrectAxSize(ax4)

 
  #plt.savefig('AS_Summary_plots.pdf', bbox_inches='tight')
  plt.show()
  plt.clf()


  xVarName = []
  for i in range(1,34):
    xVarName += [r'$A_{' + str(i) + '}$']
  indices = objAD.xWeight
  print(indices)
  dimension = len(indices)
  threshold = 1e-5

  asDict = {}
  # Add first-order indices
  for i in range(dimension):
    myKey           = '$' + xVarName[i] + '$'
    myValue         = indices[i]
    asDict[myKey]   = myValue

  # ----- Sort in decreasing order
  xKey    = np.array(asDict.keys())
  xVal    = np.array(asDict.values())
  xValAbs = np.abs(np.array(asDict.values()))
  mySort  = np.argsort(xValAbs)[::-1]
  sKey    = xKey[mySort]
  sVal    = xVal[mySort]

  # ----- Keep only the ones summing to the chosen fraction of variance
  mySum = 0.0
  majI  = len(sKey)
  for i in range(len(sKey)):
    if mySum < (1.0-threshold):
      mySum += sVal[i]*sVal[i]
    else:
      majI = i
      break

  # ----- Plot the significant indices in decreasing order

  import matplotlib as mpl
  font = {'family' : 'sans',
          'serif': 'helvet',
          'weight' : 'bold',
          'size'   : 20}
  mpl.rc('font', **font)
  mpl.rc('text', usetex=True)

  #Vanilla plot
  fig, ax1 = plt.subplots(1,1)
  ax1.scatter(range(majI), sVal[:majI], color='k', marker='o', facecolors='none')

  #ax1.set_ylim([1e-4, 1e0])
  ax1.set_xticks(range(majI))
  ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

  ax1.set_ylabel('Coefficients')
  #plt.yscale('log')

  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)

  plt.show()

  fig, ax1 = plt.subplots(1,1)
  ax1.scatter(range(majI), sVal[:majI]*sVal[:majI], color='k', marker='o', facecolors='none')

  ax1.set_ylim([1e-8, 1e0])
  ax1.set_xticks(range(majI))
  ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

  ax1.set_ylabel('Squared coefficients')
  plt.yscale('log')

  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)

  plt.show()
