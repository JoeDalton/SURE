#!/usr/bin/env python2.7

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import h5py
import math
import os
import sys
sureHome = os.environ['SURE_HOME']
sys.path.insert(0, sureHome)
import openturns as ot
from DimensionReduction.ActiveDirection import *
from SurrogateModelling.PCE import *
from Sampling.MC import *
from Utils.I_O import prt
from Numerics.Algebra import Orthonormalize
from Numerics.Algebra import VectorProjection
from Numerics.Algebra import norm2
import matplotlib as mpl
font = {'family' : 'sans',
   'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)




################################################################
# STEP 0: Important parameters
############################################################
fieldSize     = 20          # Number of "features", or field values present in the initial data
nSample       = 1000       # Number of observations, or samples, considered in the study
PCAEpsilon    = 1e-06       # Portion of the variance of the dataset ignored by the PCA
uncertainDim  = 33          # Number of uncertain input variables
nReSample     = 5000       # Number of samples of the surrogate
chosenT       = 913.5       # Temperature at which the pdf is reconstructed
UncRedEpsilon = 1e-02       # Portion of the variance ignored by the uncertainty reduction
reprSample    = 1           # Sample to reproduce on the range of temperatures



################################################################
# Some functions
############################################################
def ConstructPCE(ident, c, n, currentSamples, evals, nFeature, maxNPoly):
  thisPCE = PCE(None, ident)
  thisPCE.verbosity = False #True
  thisPCE.DefinePolyBasis(dim=c+1, quasiNorm=0.4, weights=None) 
  #thisPCE.DefineTruncatureStrat(strategy="Cleaning", maxNumPolys = maxNPoly, mostSignificant=nFeature, significanceFactor=1.0e-4)
  thisPCE.DefineTruncatureStrat(strategy="Cleaning", maxNumPolys = maxNPoly, mostSignificant=nFeature, significanceFactor=1.0e-1)
  #thisPCE.DefineEvalStrat(evalStrategy="Regression", validStrategy="KFold")
  thisPCE.DefineEvalStrat(evalStrategy="Regression", validStrategy="LOO")
  otSamples       = ot.Sample(n, c+1)
  otSamples[:,:]  = currentSamples[:,:]
  otEvals         = ot.Sample(n,1)
  for j in range(n):
    otEvals[j,0]  = evals[j]
  thisPCE.Build(otSamples, otEvals)
  return thisPCE

def LoadPCE(ID):
  thisPCE = PCE(None, ID)
  thisPCE.verbosity = False #True
  thisPCE.Load('Surrogate_PCE_' + ID + '.h5')
  return thisPCE

def SummaryPlot(data, remainder, mode, dim):
  x1 = np.linspace(data.min(), data.max(), 1000)
  fig = plt.figure(1, figsize=(7,5))
  subplots_adjust(wspace=0.2, hspace=0.35)
  ax1 = fig.add_subplot(111)
  #ax1.scatter(data, 0.5*(data + remainder), color='gray', alpha = 0.5, marker=',', s=1, label='PCE surrogate') #Seems to work for mode 2. The PCE is two times the ideal fit. strange... Idem for mode 3 but not for modes 0 and 1...
  #ax1.scatter(data, remainder, color='gray', alpha = 0.5, marker=',', s=1, label='PCE surrogate') #test
  ax1.scatter(data, (data + remainder), color='gray', alpha = 0.5, marker=',', s=1, label='PCE surrogate') 
  ax1.plot(x1,x1, color='k', label='Ideal')
  ax1.set_xlabel("Model evaluation")
  ax1.set_ylabel("Surrogate evaluation")
  ax1.set_title("Summary plot (PCE, Mode " + str (mode) + ", " + str(dim) + " dimensions")
  ax1.legend(loc='best')
  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)
  plt.savefig('testPlots/Summary_plot_Mode=' + str(mode) + '_Dim=' + str(dim) + '.pdf', bbox_inches='tight')
  plt.close()

def ModelPlot(data, surrogate, mode, abscissa):
  x = np.linspace(abscissa.min(), abscissa.max(), 1000)
  xModel = np.zeros(len(x))
  for i in range(len(x)):
    xModel[i] = surrogate.Evaluate([x[i]])[0]
  fig = plt.figure(1, figsize=(7,5))
  subplots_adjust(wspace=0.2, hspace=0.35)
  ax1 = fig.add_subplot(111)
  ax1.scatter(abscissa, data, color='gray', alpha = 0.5, marker=',', s=1, label='original data') 
  ax1.plot(x,xModel, color='k', label='Model')
  ax1.set_xlabel("Model evaluation")
  ax1.set_ylabel("Surrogate evaluation")
  ax1.set_title("Model plot (PCE, Mode " + str (mode) + ", 1 dimensions")
  ax1.legend(loc='best')
  for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
    item.set_fontsize(18)
  plt.savefig('testPlots/Model_plot_Mode=' + str(mode) + '_Dim=1.pdf', bbox_inches='tight')
  plt.close()
  print(surrogate.responseSurface)


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

#TODO: Clean the creation of xVect. this is ugly...
xVect = [np.zeros(fieldSize)] * nSample
xVect = np.array(xVect).transpose()

for i in range(fieldSize):
  T = tList[i]
  tFolderName	= '../T_' + str(T)
  #print(tFolderName)
  resultFile = tFolderName + '/result.h5'
  result      = h5py.File(resultFile, 'r')
  xEval       = np.log(result["Data/QoI"][:nSample])
  xSample     = result["Data/NormalInputVars"][:nSample,:]
  result.close()
  for j in range(nSample):
    xVect[i,j] = xEval[j]



#TODO: Here, one should center and reduce the data for better numerical behaviour !

M = xVect.T

# singular value decomposition factorises your data matrix such that:
# 
#   M = U*S*V.T     (where '*' is matrix multiplication)
# 
# * U and V are the singular matrices, containing orthogonal vectors of
#   unit length in their rows and columns respectively.
#
# * S is a diagonal matrix containing the singular values of M - these 
#   values squared divided by the number of observations will give the 
#   variance explained by each PC.
#
# * if M is considered to be an (observations, features) matrix, the PCs
#   themselves would correspond to the rows of S^(1/2)*V.T. if M is 
#   (features, observations) then the PCs would be the columns of
#   U*S^(1/2).
#
# * since U and V both contain orthonormal vectors, U*V.T is equivalent 
#   to a whitened version of M.
#
# * PCs (Modes) are already sorted by descending order 
#   of the singular values (i.e. by the
#   proportion of total variance they explain)



################################################################
# STEP 2: Computing the SVD of the data
############################################################
print("")
print("--------------------------------------------------------------")
prt('Computing the PCA of the samples', 'yellow', True)
U, s, Vt = np.linalg.svd(M, full_matrices=False)
V = Vt.T


################################################################
# STEP 3: Cutting the PCA to retrieve a certain amount of variance
############################################################
PCACutIdx = 0
threshold = (1.0-PCAEpsilon)*sum(np.square(s[:]))
for i in range(nSample):
  if sum(np.square(s[:i+1])) >= threshold:
    PCACutIdx = i
    break
prt(" Variance explained: %.6G" %(sum(np.square(s[:i+1]))/sum(np.square(s[:])) * 100) + "%", 'green', True)
prt(" " + str(PCACutIdx) + " modes used", 'green', True)


################################################################
# STEP 3bis: Sanity check 1: Plotting the "real" and reconstructed first sample to show everything is fine
############################################################
#print("")
#print("--------------------------------------------------------------")
#prt('Plotting results', 'blue', True)
nominal = M[reprSample]
reconstructed = np.zeros_like(nominal)
#for i in range(len(s)): #Exact reconstruction
for i in range(PCACutIdx): # Reduced reconstruction
  for j in range(len(reconstructed)):
    coeff = nominal.dot(V[:,i])
    reconstructed[j] += coeff * V[j,i]
#--------------------------------------------------------------------- Plot
#fig = plt.figure(figsize=(7,5))
#ax1 = fig.add_subplot(111)
##ax1.scatter(1000.0/tList, np.exp(nominal), color='k', marker='3', label=r'nominal')
#ax1.plot(1000.0/tList, np.exp(nominal), color='k', label=r'nominal')
#ax1.plot(1000.0/tList, np.exp(reconstructed), color='k', linestyle='--', label=r'reconstructed')
#ax1.set_xlabel(r'$\frac{1000}{T}$')
#ax1.set_ylabel(r'QoI')
#plt.legend(loc='best')
#for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
#  item.set_fontsize(18)
#ax1.grid()
#plt.yscale('log')
#plt.show()
#plt.close()
#fig = None


#################################################################
## STEP 3 ter: Sanity check 2: Plotting the distribution of samples projected on the PCA
#############################################################
#print("")
#print("--------------------------------------------------------------")
#prt('Plotting Results', 'blue', True)
#tIndex          = list(tList).index(chosenT)
#xReconstructed  = []
#
#for sampleIdx in range(nSample):
#  nominal = M[sampleIdx]
#  reconstructed = 0.0
#  for modeIdx in range(PCACutIdx):
#    coeff = nominal.dot(V[:,modeIdx])
#    reconstructed += coeff * V[tIndex, modeIdx]
#  xReconstructed.append(reconstructed)
#
#fig = plt.figure(figsize=(14,5))
#ax1 = fig.add_subplot(111)
#
#ax1.hist(xVect[tIndex,:], 70, color= 'dimgrey', label=r'Original', alpha = 0.5, density=True)
#ax1.hist(xReconstructed, 70, label=r'PCA Reconstruction with ' + str(PCACutIdx +1) + ' modes', alpha = 0.5, density=True)
#
#ax1.set_xlabel(r'log(IDT (s))')
#ax1.set_ylabel(r'PDF')
#
##plt.xscale('log')
#plt.title('KDE of the original and reconstructed pdfs at T = 913.5K')
#plt.legend(loc='best')
#for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
#  item.set_fontsize(18)
#plt.show()
##quit()

################################################################
# STEP 4: Building the (new) reduced datasets (projections of the data on the modes)
############################################################
print("")
print("--------------------------------------------------------------")
prt('Building reduced datasets', 'yellow', True)
xDataset = []
for i in range(PCACutIdx):
  dataset = np.zeros(nSample)
  for j in range(nSample):
    dataset[j] = M[j].dot(V[:,i])
  xDataset.append(dataset)
#TODO: Center and reduce the datasets to improve numerical behaviour. Warning, reconstruction must also be prepared...


################################################################
# STEP 5: Building PCE surrogates for each dataset
############################################################
xPCE      = []
xxAD      = []
for i in [2]: #range(PCACutIdx): #DEBUG
  xAD = []
  xxRedSample = []
  print("")
  print("--------------------------------------------------------------")
  prt('Building a PCE surrogate for mode ' + str(i), 'yellow', True)

  #Step 0. Initialization
  data          = np.zeros(nSample)
  data[:]       = xDataset[i][:]
  counter       = 0
  xRemInput     = np.zeros_like(xSample)
  xRemInput[:]  = xSample[:]
  
  while counter < 3: #uncertainDim:
    prt(' Finding active direction ' + str(counter) +'...', 'green', True)
    
    #print('test')
    #if counter>0:
    #  for sampleIdx in range(nSample):
    #    #print(xAD[-1].xWeight[:].dot(xRemInput[sampleIdx,:]) * xAD[-1].xWeight[:])
    #    #print(xAD[-1].xWeight[:].dot(xRemInput[sampleIdx,:]))
    #    print(xRemInput[sampleIdx,:].dot(xAD[-1].xWeight) / norm2(xRemInput[sampleIdx,:]))
    #  break


    # Step 1. (Finding the AD of the remainder)
    ID = 'Mode_' + str(i) + '_RedDim_' + str(counter)
    AD = ActiveDirection(ID, xRemInput, data, parent=None, verbosity=False) #True)
    xAD.append(AD)

    # DEBUG
    if counter >= 1:
      print('Cos of the two previously calculated directions = ' + str(xAD[counter-1].xWeight.dot(xAD[counter].xWeight)/ (norm2(xAD[counter-1].xWeight) * norm2(xAD[counter].xWeight) )))
      #quit()

    # Step 1.1.(Orthonormalizing the new AD with regard to the basis)
    prt(' Orthonormalizing the basis...', 'green', True)
    xTempWeight = []
    for k in range(counter+1):
      xTempWeight.append(xAD[k].xWeight)
    Orthonormalize(xTempWeight)
    #xTempWeight = Orthonormalize(xTempWeight)
    xAD[-1].xWeight = xTempWeight[-1] # TODO: Maybe do this in a cleaner way
   
    # DEBUG
    if counter >= 1:
      print('Cos of the two last orthogonal directions = ' + str(xAD[counter-1].xWeight.dot(xAD[counter].xWeight)/ (norm2(xAD[counter-1].xWeight) * norm2(xAD[counter].xWeight) )))
      #quit()

    #print("----------------------------------------------------------------------")
    newXRemInput = np.zeros_like(xRemInput)
    for sampleIdx in range(nSample):
      #xRemInput[sampleIdx,:] = xRemInput[sampleIdx,:] - VectorProjection(xRemInput[sampleIdx,:], xAD[-1].xWeight)
      newXRemInput[sampleIdx,:] = xRemInput[sampleIdx,:] - VectorProjection(xRemInput[sampleIdx,:], xAD[-1].xWeight)
      print(VectorProjection(xRemInput[sampleIdx,:], xAD[-1].xWeight).dot(newXRemInput[sampleIdx,:]))
    
    xRemInput[:] = newXRemInput[:]

    ## Step 1.2.(Preparing the new reduced samples)
    #xRedSample = AD.ComputeReducedFromTotal(xSample)
    #xxRedSample.append(xRedSample)
    #xCurrentSample = np.zeros((nSample, counter+1))
    #for j in range(nSample):
    #  for k in range(counter+1):
    #    xCurrentSample[j,k] = xxRedSample[k][j]

    ## Step 2. (Computing the PCE)
    #prt(' Computing the PCE over the ' + str(counter+1) + ' first directions...', 'green', True)
    #maxNPoly = 500
    ##nFeature = int(3 * (counter+1) * (floor(math.log(counter+2)) + 1))
    #nFeature = 50 
    #myPCE = ConstructPCE(ID, counter, nSample, xCurrentSample, xDataset[i], nFeature, maxNPoly)   


    ##Step 3. (Extracting the remainder)
    #xRemain = np.zeros(nSample)
    #for j in range(nSample):
    #  xRemain[j] = data[j] - myPCE.Evaluate(xCurrentSample[j,:])

    ##Step 3.1. (Computing the remainding variance)
    #testRatio = np.var(xRemain)/np.var(xDataset[i])
    #print("Proportion of the variance discarded: %.3G" %(testRatio * 100) + '%')
   
    ##Step 3.2.(Plotting a summary plot)
    #SummaryPlot(data, data + xRemain, i, counter+1)
    ##Step 3.3 (Plotting a model plot when there is 1 dim)
    #if counter ==0:
    #  ModelPlot(data, myPCE, i, xCurrentSample[:,0])
    #  #quit()
  
    ##Step 4. Looping or not looping... 
    counter += 1
    #if testRatio < UncRedEpsilon:
    #  break


  # TODO: save all the AD corresponding to the PCEs !
  #myPCE.Save()
  #myPCE.Dump()
 
  ##Loading the PCEs
  #myPCE = LoadPCE(ID)

  #xPCE.append(myPCE)

quit() #DEBUG









#################################################################
## STEP 6: Plotting again the "real" and reduced + (reduced +) surrogated first sample to show everything is fine
#############################################################
#nominal = M[reprSample]
#reconstructed2 = np.zeros_like(nominal)
#for i in range(PCACutIdx): # Reduced reconstruction
#  for j in range(len(reconstructed)):
#    # For each mode, the coeff must be reconstructed from its surrogate
#    coeff = xPCE[i].Evaluate(xSample[reprSample])
#    reconstructed2[j] += coeff * V[j,i]
##--------------------------------------------------------------------- Plot
#print("")
#print("--------------------------------------------------------------")
#prt('Plotting results', 'blue', True)
#fig = plt.figure(figsize=(14,10))
#ax1 = fig.add_subplot(111)
##ax1.scatter(1000.0/tList, np.exp(nominal), color='k', marker='3', label=r'nominal')
#ax1.plot(1000.0/tList, np.exp(nominal), color='k', label=r'nominal')
#ax1.plot(1000.0/tList, np.exp(reconstructed), color='k', linestyle='--', label=r'PCA reconstruction')
#ax1.plot(1000.0/tList, np.exp(reconstructed2), color='k', linestyle=':', label=r'PCA+AD+spline recon.')
#ax1.set_xlabel(r'$\frac{1000}{T}$')
#ax1.set_ylabel(r'QoI')
#plt.legend(loc='best')
#for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
#  item.set_fontsize(18)
#ax1.grid()
#plt.yscale('log')
##plt.savefig('PCA+PCEReconstruction.pdf', bbox_inches='tight')
#plt.show()
#plt.close()


################################################################
# STEP 7: Resampling the reduced + (reduced +) surrogates at the chosen temperature
############################################################
print("")
print("--------------------------------------------------------------")
prt('Resampling for T = 913.5K', 'yellow', True)
distribution  = ot.ComposedDistribution([ot.Normal(0.0,1.0)] * uncertainDim)
MCDrawer      = MC(None, distribution, nReSample, firstSample=0, verbosity=False, fixSeed=False)
xReSample     = MCDrawer.Draw()
#QMCDrawer     = QMC(None, distribution, nReSample, firstSample=0, verbosity=False, sequence='Sobol')
#xReSample     = QMCDrawer.Draw()
xReEval       = np.zeros(nReSample)
tIndex        = list(tList).index(chosenT)

for modeIdx in range(PCACutIdx):
  for sampleIdx in range(nReSample):
    coeff        = xPCE[modeIdx].Evaluate(xReSample[sampleIdx])
    xReEval[sampleIdx] += coeff * V[tIndex, modeIdx]


#################################################################
## STEP 8: Computing the Gaussian KDE estimations of the initial and reconstructed pdfs
#############################################################
#print("")
#print("--------------------------------------------------------------")
#prt('Building KDEs of the original and reconstructed samples', 'yellow', True)
#nDrawingPoint   = 1000
#xDrawingPoint   = np.linspace(xVect[tIndex,:].min(), xVect[tIndex,:].max(), nDrawingPoint)
#
#def normPDF(x):
#  return (1.0/(math.sqrt(2.0*math.pi))) * math.exp(-0.5*x*x)
#
#def GaussianKernelPDF(x, arr, h):
#  temp = 0
#  for i in range(len(arr)):
#    temp += normPDF((x-arr[i])/h)
#  temp /= h
#  temp /= len(arr)
#  return temp
#
#xOrigPDF    = np.zeros(nDrawingPoint)
#xRecoPDF    = np.zeros(nDrawingPoint)
#
## Finding optimum bandwith for the estimators
#from sklearn.neighbors import KernelDensity
#from sklearn.model_selection import GridSearchCV
#params = {'bandwidth': np.logspace(-1, 1, 20)}
##grid = GridSearchCV(KernelDensity(), params, cv=5, iid=False, n_jobs=-1) #5-fold cross-validation - n_jobs=-1 => All cores
#grid = GridSearchCV(KernelDensity(), params, cv=None, iid=False, n_jobs=-1) # Default 3-fold cross validation
#
#prt('Optimizing bandwith for original samples', 'green', True)
#grid.fit(xVect[tIndex,:].reshape(-1,1))
#hOrig = grid.best_estimator_.bandwidth
#
#prt('Optimizing bandwith for reconstructed samples', 'green', True)
#grid.fit(xReEval.reshape(-1,1))
#hReco = grid.best_estimator_.bandwidth
#
## Computing the KDE
#prt('Computing the KDE estimators', 'green', True)
#for i in range(nDrawingPoint): 
#  xOrigPDF[i]   = GaussianKernelPDF(xDrawingPoint[i], xVect[tIndex,:],  hOrig)
#  xRecoPDF[i]   = GaussianKernelPDF(xDrawingPoint[i], xReEval,          hReco)


################################################################
# STEP 9: Plotting the estimated pdfs
############################################################
print("")
print("--------------------------------------------------------------")
prt('Plotting Results', 'blue', True)
fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(111)


ax1.hist(xVect[tIndex,:], 70, label=r'Original', alpha = 0.5)
ax1.hist(xReEval, 70, label=r'Reconstructed', alpha = 0.5)
#ax1.plot(xDrawingPoint, xOrigPDF, color='k', linestyle='-',label=r'Original')
#ax1.plot(xDrawingPoint, xRecoPDF, color='k', linestyle=':',label=r'Reconstructed')

ax1.set_xlabel(r'log(IDT (s))')
ax1.set_ylabel(r'PDF')

#plt.xscale('log')
plt.title('KDE of the original and reconstructed pdfs at T = 913.5K')
plt.legend(loc='best')
for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(18)
plt.show()
#plt.savefig('KDE_PCE.pdf', bbox_inches='tight')
