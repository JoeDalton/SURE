import numpy as np
import sys
from DimensionReduction  import *
from Utils.I_O import prt
from Utils.I_O import ProgBar
from Numerics.Algebra import norminf
from Numerics.Algebra import norm1
from Numerics.Algebra import norm2
from Numerics.Statistics import ComputeNormalizedDistanceFromSamples as ComputeDist
from Preprocessing import CenterScaling

def FindOptimalPCA(threshold, data, PCAID='myPCA', CSID="myCS", parent=None, verbosity=-1, norm='norminf'):

  # Should the error really be normalized ?
  

  myNModes        = {}
  myErrors        = {}
  nFeature        = data.shape[0]
  nSample         = data.shape[1]
  currentOptNMode = nFeature
  
  methodList = ['VAST', 'AUTO', 'RANGE', 'PARETO', 'LEVEL', 'MAX']
  #methodList = ['PARETO']
  
  if verbosity >= 0:
    myProgBar = ProgBar('Looping over scaling methods', len(methodList))


  ################################################################
  # STEP 0: Looping over the scaling methods
  ############################################################
  for scaleMethod in methodList:
  
    ################################################################
    # STEP 1: Preprocessing samples
    ############################################################
    #myCS = CenterScaling(CSID, data, scaleMethod, isCGlobal=False, isSGlobal=False)
    myCS = CenterScaling(CSID, data, scaleMethod, isCGlobal=True, isSGlobal=True)
    # TODO: Also check global scaling ?
    pro = myCS.ComputePro(data)
  
    ################################################################
    # STEP 2: Looping over the number of modes
    ############################################################
    prt("", 'None', verbosity)
    prt("--------------------------------------------------------------", 'None', verbosity-1)
    prt('Building the PCAs using the ' + scaleMethod + ' scaling method', 'yellow', verbosity-1)

    nMode = 1
    while nMode <= currentOptNMode:
      ################################################################
      # STEP 3: Building the PCA of the data
      ############################################################
      myPCA = PCA(PCAID, pro, nMode=nMode, parent=None, verbosity=verbosity-1)
      
      ################################################################
      # STEP 4: Computing the error
      ############################################################
      prt("--------------------------------------------------------------", 'None', verbosity-1)
      prt('Computing the error', 'None', verbosity-1)
     
      redMat          = myPCA.ComputeReducedFromTotal(pro.T)
      reconstructed   = myCS.ComputeOrig(myPCA.ComputeTotalFromReduced(redMat).T).T
      testVect = []
      for featIdx in range(nFeature): 
        testVect.append(ComputeDist(data[featIdx,:], reconstructed[:,featIdx]))
      if norm == 'norminf':
        error = norminf(testVect)
      elif norm == 'norm2':
        error = norm2(testVect)
      elif norm == 'norm1':
        error = norm1(testVect)
      else:
        prt('', 'red', True)
        prt('DimensionReduction/PCA.py/FindOptimalPCA', 'red', True)
        prt('Error: Please provide a valid "norm" argument. Valid arguments are "norminf", "norm2" and "norm1".', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit() 

      myPCA = None
      if error < threshold:
        prt('Current error = ' + str(error), 'None', verbosity-1)
        if nMode < currentOptNMode:
          currentOptNMode = nMode
        break
      else:
        prt('Current error = ' + str(error), 'None', verbosity-1)
        nMode += 1
    print("") #DEBUG
    print(scaleMethod + ': nMode = ' + str(nMode) + ', error = ' + str(error)) #DEBUG
    if verbosity >= 0:
      myProgBar.Update()
    myNModes[scaleMethod] = nMode
    myErrors[scaleMethod] = error # Will be used to sort out ex-aequo winners


  
  ################################################################
  # STEP 5: Finding the optimal scaling
  ############################################################
  # Sorting out ex-aequo winners on the num. of modes
  xTempWinner = []
  optNMode = nFeature
  for scaleMethod in methodList:
    if myNModes[scaleMethod] == optNMode:
      xTempWinner.append(scaleMethod)
    elif myNModes[scaleMethod] < optNMode:
      xTempWinner = []
      xTempWinner.append(scaleMethod)
      optNMode = myNModes[scaleMethod]
    else:
      # Method goes to garbage
      pass

  # Getting the one with the lowest error
  minError = 1.0e60
  #optMethod = ''
  for scaleMethod in xTempWinner:
    if myErrors[scaleMethod] < minError:
      minError = myErrors[scaleMethod]
      optMethod = scaleMethod


  #optMethod = min(myNModes, key=myNModes.get) # Must be updates to allow ex-aequo winners sorted by their error
  optNMode  = myNModes[optMethod]
  
  prt('', 'None', verbosity-1)
  prt('------------------------------', 'None', verbosity-1)
  prt('The winner is ' + optMethod + ' with ' + str(optNMode) + ' modes.', 'green', verbosity-1)

  optCS     = CenterScaling(CSID, data, optMethod, isCGlobal=False, isSGlobal=False)
  pro       = optCS.ComputePro(data)
  optPCA    = PCA(PCAID, pro, nMode=optNMode, parent=parent, verbosity=verbosity)
 
  if verbosity >= 0:
    myProgBar.Terminate()
  return optCS, optPCA




class PCA(DimensionReduction):
  """
  PCA is computed with the SVD method
  Singular Value Decomposition factorises your data matrix such that:
   
     M = U*S*V.T     (where '*' is matrix multiplication)
   
   * U and V are the singular matrices, containing orthogonal vectors of
     unit length in their rows and columns respectively.
  
   * S is a diagonal matrix containing the singular values of M - these 
     values squared divided by the number of observations will give the 
     variance explained by each PC.
  
   * if M is considered to be an (observations, features) matrix, the PCs
     themselves would correspond to the rows of S^(1/2)*V.T. if M is 
     (features, observations) then the PCs would be the columns of
     U*S^(1/2).
  
   * since U and V both contain orthonormal vectors, U*V.T is equivalent 
     to a whitened version of M.
  
   * PCs (Modes) are already sorted by descending order 
     of the singular values (i.e. by the
     proportion of total variance they explain)
  """
  ######################
  #     Properties     #
  ######################
  epsilon            = 0.0
  xMode              = None
  xSing              = None

  ###################
  #     Methods     #
  ###################
  def __init__(self, ID, data, parent=None, verbosity=-1, **kwargs):
    # data is of shape (nFeature, nSample)
    # self.xMode is of shape (nMode, nFeature)
    # Initializing
    self.ID                = ID 
    self.parent            = parent
    self.verbosity         = verbosity
    self.savedItems        = {}
    self.totalDim          = data.shape[0]
    self.objectType        = 'PCA'
    if "epsilon" in kwargs and "nMode" not in kwargs:
      self.epsilon         = kwargs["epsilon"]
    elif "nMode" in kwargs and "epsilon" not in kwargs:
      self.reducedDim      = kwargs["nMode"]
      self.epsilon         = 0.0
    else:
      prt('', 'red', True)
      prt('DimensionReduction/PCA.py/PCA/__init__', 'red', True)
      prt('Error: Either "epsilon" or "nMode" must be provided as kwargs. They can not be provided together.', 'red', True)
      prt('Exiting', 'red', True)
      sys.exit() 

    #Computing the PCA
    prt("", 'green', self.verbosity)
    prt("Computing PCA " + ID + "....", 'green', self.verbosity)
    self.xMode                = None
    self.xWeight              = None
    U, self.xSing, self.xMode = np.linalg.svd(data.T, full_matrices=False)

    if "epsilon" in kwargs:
      nMode                     = 0
      threshold                 = (1.0-self.epsilon) * sum(np.square(self.xSing[:]))
      nSample                   = data.shape[1]
      for i in range(nSample):
        if sum(np.square(self.xSing[:i+1])) >= threshold:
          nMode = i
          break
        if i == nSample-1:
          nMode = i
      self.reducedDim = nMode

    prt("PCA " + ID + " computed.", 'green', self.verbosity)
    prt("Variance explained: %.6G" %(sum(np.square(self.xSing[:self.reducedDim+1]))/sum(np.square(self.xSing[:])) * 100) + "%", 'green', self.verbosity)
    prt("" + str(self.reducedDim) + " modes used", 'green', self.verbosity)


  def Save(self):
    self.savedItems['ID']         = self.ID
    self.savedItems['Verbosity']  = self.verbosity
    self.savedItems['TotalDim']   = self.totalDim
    self.savedItems['ReducedDim'] = self.reducedDim
    self.savedItems['ObjectType'] = "Reduction_" + self.objectType
    self.savedItems['xSing']      = self.xSing
    self.savedItems['xMode']      = self.xMode


  def ComputeReducedFromTotal(self, totalVar):
    if totalVar.ndim == 1:
      if self.totalDim == len(totalVar):
        redVar = np.zeros(self.reducedDim)
        for i in range(self.reducedDim):
          redVar[i] = totalVar.dot(self.xMode[i,:])
        return redVar
      else:
        prt('', 'red', True)
        prt('DimensionReduction/PCA.py/PCA/ComputeReducedFromTotal', 'red', True)
        prt('Error: The dimension of totalVar is not equal to the total dimension of the reduction method', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit() 

    elif totalVar.ndim == 2:
      nSample   = totalVar.shape[0]
      sampleDim = totalVar.shape[1]
      if self.totalDim == sampleDim:
        xRedVar = np.zeros((nSample, self.reducedDim))
        for j in range(nSample):
          for modeIdx in range(self.reducedDim):
            xRedVar[j,modeIdx] = totalVar[j].dot(self.xMode[modeIdx,:])
        return xRedVar
      else:
        prt('', 'red', True)
        prt('DimensionReduction/PCA.py/PCA/ComputeReducedFromTotal', 'red', True)
        prt('Error: The dimension of totalVar[i] is not equal to the total dimension of the reduction method', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit() 

    else:
      prt('', 'red', True)
      prt('DimensionReduction/PCA.py/PCA/ComputeReducedFromTotal', 'red', True)
      prt('Error: Incompatible shape of totalVar', 'red', True)
      prt('Exiting', 'red', True)
      sys.exit() 
  

  def ComputeTotalFromReduced(self, reducedVar):

    if reducedVar.ndim == 1:
      if self.reducedDim == len(reducedVar):
        totVar = np.zeros(self.totalDim)
        for modeIdx in range(self.reducedDim):
          #for featIdx in range(self.totalDim): #TODO: vectorise for better performance ?
          #  totVar[featIdx] += reducedVar[modeIdx] * self.xMode[modeIdx,featIdx]
          totVar[:] += reducedVar[modeIdx] * self.xMode[modeIdx,:]
        return totVar
      else:
        prt('', 'red', True)
        prt('DimensionReduction/PCA.py/PCA/ComputeTotalFromReduced', 'red', True)
        prt('Error: The dimension of reducedVar is not equal to the reduced dimension of the reduction method', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit() 

    elif reducedVar.ndim == 2:
      nSample   = reducedVar.shape[0]
      sampleDim = reducedVar.shape[1]
      if self.reducedDim == sampleDim:
        xTotVar = np.zeros((nSample, self.totalDim))
        for sampleIdx in range(nSample):
          for modeIdx in range(self.reducedDim):
            for featIdx in range(self.totalDim):
              xTotVar[sampleIdx,featIdx] += reducedVar[sampleIdx,modeIdx] * self.xMode[modeIdx,featIdx]
        return xTotVar
      else:
        prt('', 'red', True)
        prt('DimensionReduction/PCA.py/PCA/ComputeTotalFromReduced', 'red', True)
        prt('Error: The dimension of reducedVar[i] is not equal to the reduced dimension of the reduction method', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit() 

    else:
      prt('', 'red', True)
      prt('DimensionReduction/PCA.py/PCA/ComputetotalFromReduced', 'red', True)
      prt('Error: Incompatible shape of totalVar', 'red', True)
      prt('Exiting', 'red', True)
      sys.exit() 
