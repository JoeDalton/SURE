#!/usr/bin/env python2.7


import numpy as np
import matplotlib.pyplot as plt
import Link
import openturns as ot
from Sampling.MC   import *
from Numerics.Algebra import Orthonormalize
from Numerics.Algebra import VectorProjection
from Numerics.Algebra import norm2
from DimensionReduction.ActiveDirection import *
from SurrogateModelling.PCE import *
import matplotlib as mpl
font = {'family' : 'sans',
   'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)



#------------------------------------------------------------ Meta-parameters
nSample      = 500
firstSample  = 0
nDim         = 3
distribution = ot.ComposedDistribution([ot.Uniform(0.0, 1.0)] * nDim)
myThreshold  = 1.0e-4






#------------------------------------------------------------ Preparation of samples
MCDrawer   = MC(None, distribution, nSample, firstSample=firstSample, verbosity=False, fixSeed=False)
s1 = MCDrawer.Draw()
#x1 = np.zeros(nSample)
#y1 = np.zeros(nSample)
#x1 = s1[:, 0]
#y1 = s1[:, 1]







#------------------------------------------------------------ Miscellaneous plots
def ReductionPlot(data, abscissa):
  x = np.linspace(abscissa.min(), abscissa.max(), 1000)
  fig = plt.figure(1, figsize=(7,5))
  ax1 = fig.add_subplot(111)
  ax1.scatter(abscissa, data, color='gray', alpha = 0.5, marker=',', s=1, label='original data') 
  ax1.set_xlabel("reduced variable")
  ax1.set_ylabel("data")
  ax1.set_title("Reduction plot")
  ax1.legend(loc='best')
  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)
  plt.show()
  plt.close()

def SpacePlot(xHyperplanPoint, xDirectionPoint):
  #from matplotlib import animation
  from mpl_toolkits.mplot3d import Axes3D
  fig = plt.figure()
  ax = Axes3D(fig)
  
  def init():
      ax.scatter(xHyperplanPoint[:,0], xHyperplanPoint[:,1], xHyperplanPoint[:,2], marker='3', s=20, c="k", alpha=0.6, label='Projected on hyperplan')
      ax.scatter(xDirectionPoint[:,0], xDirectionPoint[:,1], xDirectionPoint[:,2], marker='o', s=20, c="k", alpha=0.6, label='Projected on direction')
      return fig,
  plt.legend(loc='best')
  init()
  plt.show()
  plt.close()
  fig = None


def SummaryPlot(data, remainder, dim):
  x1 = np.linspace(data.min(), data.max(), 1000)
  fig = plt.figure(1, figsize=(7,5))
  ax1 = fig.add_subplot(111)
  ax1.scatter(data, (data + remainder), color='gray', alpha = 0.5, marker=',', s=1, label='PCE surrogate') 
  ax1.plot(x1,x1, color='k', label='Ideal')
  ax1.set_xlabel("Model evaluation")
  ax1.set_ylabel("Surrogate evaluation")
  ax1.set_title("Summary plot (PCE, " + str(dim) + " dimensions")
  ax1.legend(loc='best')
  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)
  #plt.savefig('testPlots/Summary_plot_Dim=' + str(dim) + '.pdf', bbox_inches='tight')
  plt.show()
  plt.close()


def ConstructPCE(ident, c, n, currentSamples, evals, nFeature, maxNPoly):
  thisPCE = PCE(None, ident)
  thisPCE.verbosity = False #True
  thisPCE.DefinePolyBasis(dim=c+1, quasiNorm=0.4, weights=None) 
  thisPCE.DefineTruncatureStrat(strategy="Cleaning", maxNumPolys = maxNPoly, mostSignificant=nFeature, significanceFactor=1.0e-4)
  #thisPCE.DefineEvalStrat(evalStrategy="Regression", validStrategy="KFold")
  thisPCE.DefineEvalStrat(evalStrategy="Regression", validStrategy="None")
  otSamples       = ot.Sample(n, c+1)
  otSamples[:,:]  = currentSamples[:,:]
  otEvals         = ot.Sample(n,1)
  for j in range(n):
    otEvals[j,0]  = evals[j]
  thisPCE.Build(otSamples, otEvals)
  return thisPCE





#------------------------------------------------------------ Defining the model and evaluating samples
def MyModel3(x,y,z):
  result = 0.0
  #result += x**2 + 2*y**2 + 3*z**2 #+ x*y + 0.5*x*z + 4*y*z
  result += 3*(x+y)**2 + 0.5*(y+z)**2 #+ 3*z**2 #+ x*y + 0.5*x*z + 4*y*z
  return result

def MyModel2(x,y):
  result = 0.0
  result += x**2 + 2*y**2 #+ x*y
  return result


data = np.zeros(nSample)
for sampleIdx in range(nSample):
  data[sampleIdx] = MyModel3(s1[sampleIdx,0],s1[sampleIdx,1],s1[sampleIdx,2])
  #data[sampleIdx] = MyModel2(s1[sampleIdx,0],s1[sampleIdx,1])





#------------------------------------------------------------ Finding active directions
counter         = 0
xRemInput       = np.zeros((nSample,nDim))
xRemInput[:,:]  = s1[:,:]
xAD             = []
remData         = np.zeros_like(data)
remData[:]    = data[:]
xxRedVar        = []

while counter < nDim:
  prt(' Finding active direction ' + str(counter) +'...', 'None', True)
  ###################################################################################### Step 1. Finding the current AD
  ID = 'AD_' + str(counter)
  AD = ActiveDirection(ID, xRemInput, remData, parent=None, verbosity=False)
  xAD.append(AD)




  ###################################################################################### Step 2. Orthonormalizing the basis
  prt(' Orthonormalizing the basis...', 'None', True)
  xTempWeight = []
  for k in range(counter+1):
    xTempWeight.append(xAD[k].xWeight)
  xTempWeight = Orthonormalize(xTempWeight)
  #xTempWeight = Orthonormalize(xTempWeight)
  xAD[-1].xWeight = xTempWeight[-1] # TODO: Maybe do this in a cleaner way
 
  xRedVar = xAD[-1].ComputeReducedFromTotal(xRemInput)
  xxRedVar.append(xRedVar)

  #----------------------- Diagnostics
  print("Computed direction:")
  print(xAD[-1].xWeight)
  ReductionPlot(remData, xRedVar)
  #if counter >= 1:
  #  print('Cos of the two last directions:')
  #  print(xAD[-1].xWeight.dot(xAD[counter-1].xWeight))
  # ----------------------



  ###################################################################################### Step 3. Computing a surrogate and the remainder
  prt(' Computing the PCE over the ' + str(counter+1) + ' first directions...', 'green', True)
  xCurrentSample = np.zeros((nSample, counter+1)) #TODO: Check if this intermediary is really necessary ?
  for sampleIdx in range(nSample):
    for dim in range(counter+1):
      xCurrentSample[sampleIdx,dim] = xxRedVar[dim][sampleIdx]

  myPCE = ConstructPCE(ID, counter, nSample, xCurrentSample, data, 50, 500)   

  for sampleIdx in range(nSample):
    remData[sampleIdx] = remData[sampleIdx] - myPCE.Evaluate(xCurrentSample[sampleIdx,:])




  #----------------------- Diagnostics
  SummaryPlot(data, remData, counter+1)
  #if counter ==0:
  #  ModelPlot(data, myPCE, i, xCurrentSample[:,0])
  #-----------------------





  ###################################################################################### Step 4. To loop or not to loop
  testRatio = np.var(remData)/np.var(data)
  print("Proportion of the variance discarded: %.3G" %(testRatio * 100) + '%')
  if testRatio < myThreshold:
    break



  ###################################################################################### Step 5. Computing the projection of the inputs on the hyperplan for the next ite
  newXRemInput = np.zeros_like(xRemInput)
  dirProj      = np.zeros_like(xRemInput) #DEBUG
  #print("projected on direction")
  for sampleIdx in range(nSample):
    newXRemInput[sampleIdx,:] = xRemInput[sampleIdx,:] - VectorProjection(xRemInput[sampleIdx,:], xAD[-1].xWeight)
    dirProj[sampleIdx,:] = VectorProjection(xRemInput[sampleIdx,:], xAD[-1].xWeight) 
  
  #----------------------- Diagnostics
  if counter > -1:
    SpacePlot(newXRemInput, dirProj)
  #-----------------------
  

  xRemInput[:,:] = newXRemInput[:,:]
  

  counter += 1





