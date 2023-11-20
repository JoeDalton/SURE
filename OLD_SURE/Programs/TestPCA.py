#!/usr/bin/env python2.7

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import h5py
import Link
from DimensionReduction.Preprocessing import *
from DimensionReduction.PCA import *
from Utils.I_O import prt
import matplotlib as mpl #TODO: move this in a special module (Utils.I_O ?). It works and would purify the programs!
font = {'family' : 'sans',
   'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)


from Numerics.Statistics import *

################################################################
# STEP 0: Important parameters
############################################################
fieldSize     = 20          # Number of "features", or field values present in the initial data
nSample       = 65000       # Number of observations, or samples, considered in the study
PCAEpsilon    = 1.0e-2 #5e-04       # Portion of the variance of the dataset ignored by the PCA
uncertainDim  = 33          # Number of uncertain input variables
nReSample     = 65000       # Number of samples of the surrogate
chosenT       = 913.5       # Temperature at which the pdf is reconstructed
UncRedEpsilon = 1e-02       # Portion of the variance ignored by the uncertainty reduction
reprSample    = 1           # Sample to reproduce on the range of temperatures

#threshold = 5.0e-6
#threshold = 1.0e-5 # Donne OK en VAST, mais avec 7 modes (norminf), 8 modes LEVEL en norm2
#threshold = 5.0e-5 # 
threshold = 2.3e-5 # 


################################################################
# STEP 1: Loading all the samples
############################################################
print("")
print("--------------------------------------------------------------")
prt('Loading samples', 'yellow', True)
precursorList = np.linspace(0.95, 1.2, fieldSize)
tList = 1000/precursorList
tList[:] = np.around(tList[:], 1)
indices = np.linspace(1, 33 , 33)

xVect = np.zeros((fieldSize, nSample))
for i in range(fieldSize):
  T = tList[i]
  resultFile = 'testData/result_' +str(T) + '.h5'
  result      = h5py.File(resultFile, 'r')
  xEval       = np.log(result["Data/QoI"][:nSample])
  xSample     = result["Data/NormalInputVars"][:nSample,:]
  result.close()
  for j in range(nSample):
    xVect[i,j] = xEval[j]


################################################################
# STEP 2: Building the optimal couple of scaling and PCA to reach the target threshold
############################################################
myCS, myPCA = FindOptimalPCA(threshold, xVect, verbosity = False, norm='norminf')  
  
  
################################################################
# STEP 3: Plotting the "real" and reconstructed first sample
############################################################
##print("")
##print("--------------------------------------------------------------")
#prt('Plotting first sample', 'None', True)
#initTot  = xVect.T[reprSample]
#reduced  = myPCA.ComputeReducedFromTotal(myCS.ComputePro(initTot))
#reconTot = myCS.ComputeOrig(myPCA.ComputeTotalFromReduced(reduced))

##--------------------------------------------------------------------- Plot
#fig = plt.figure(figsize=(7,5))
#ax1 = fig.add_subplot(111)
#ax1.scatter(1000.0/tList, np.exp(initTot), color='k', marker='3', label=r'nominal')
#ax1.plot(1000.0/tList, np.exp(reconTot), color='k', linestyle='--', label=r'reconstructed')
#ax1.set_xlabel(r'$\frac{1000}{T}$')
#ax1.set_ylabel(r'QoI')
#plt.legend(loc='best')
#plt.title(scaleMethod + ', ' +str(myPCA.reducedDim) + ' modes')
#for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
#  item.set_fontsize(18)
#ax1.grid()
#plt.yscale('log')
##plt.show()
#plt.savefig('results/Sample_'+str(reprSample)+'_reconstruction_'+scaleMethod+'.pdf', bbox_inches='tight')
#plt.close()
#fig = None




#################################################################
## STEP 4: Plotting the distribution of samples projected on the PCA at the chosen T
#############################################################
print("")
print("--------------------------------------------------------------")

pro             = myCS.ComputePro(xVect)
redMat          = myPCA.ComputeReducedFromTotal(pro.T)
reconstructed   = myCS.ComputeOrig(myPCA.ComputeTotalFromReduced(redMat).T).T
scaleMethod     = myCS.myScaling.method

for tIndex in range(len(tList)):
  prt('Plotting The distribution for T = ' + str(tList[tIndex]) + 'K', 'None', True)
  
  D = ComputeDistanceFromSamples(xVect[tIndex,:], reconstructed[:,tIndex])
  ND = ComputeNormalizedDistanceFromSamples(xVect[tIndex,:], reconstructed[:,tIndex])
  
  
  fig = plt.figure(figsize=(7,5))
  ax1 = fig.add_subplot(111)
  
  ax1.hist(xVect[tIndex,:], 70, color= 'dimgrey', label=r'Original', alpha = 0.5, density=True)
  ax1.hist(reconstructed[:,tIndex], 70, label=r'PCA Reconstruction with ' + str(myPCA.reducedDim) + ' modes', alpha = 0.5, density=True)
  
  ax1.set_xlabel(r'log(IDT (s))')
  ax1.set_ylabel(r'PDF')
  
  plt.title(scaleMethod + ', ' + str(myPCA.reducedDim) + ' modes, D = ' + str(D) + ', ND = ' + str(ND))
  plt.legend(loc='best')
  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)
  #plt.show()
  plt.savefig('results/T_'+str(tList[tIndex])+'_reconstruction_'+scaleMethod+'.pdf', bbox_inches='tight')
  plt.close()
  fig = None



  
