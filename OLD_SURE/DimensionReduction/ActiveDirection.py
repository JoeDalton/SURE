import numpy as np
import sys
from sklearn.linear_model import LinearRegression
from DimensionReduction  import *
from Utils.I_O import prt, DictFromH5, ErrorMsg
from Numerics.Algebra import norm2

fileName = 'DimensionReduction/ActiveDirection.py'

class ActiveDirection(DimensionReduction):
  
  ######################
  #     Properties     #
  ######################
  xWeight            = 0
  determinationCoeff = 0.0

  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    self.savedItems         = {}
    self.reducedDim         = 1
    self.determinationCoeff = 0.0
    self.objectType         = 'Reduction_Preprocessing_ActiveDirection'
    willExit                = True
    if "block" in kwargs and "experiment" in kwargs:
      self.parent = kwargs["experiment"]  # parent must be an experiment object
      self.totalDim = self.parent.distribution.getDimension()
      self.xWeight = np.zeros(self.totalDim)
      self._InitWithBlock(kwargs["block"]) # kwargs["block"] must be a block of "InputHandler"
      willExit = False
    elif "block" in kwargs:
      self._InitWithBlock(kwargs["block"]) # kwargs["block"] must be a block of "InputHandler"
      willExit          = False
    elif "empty" in kwargs and "verbosity" in kwargs and "totalDim" in kwargs:
      self.verbosity    = kwargs["verbosity"]
      willExit          = False
    else:
      try:
        self.ID           = kwargs["ID"]
        if "experiment" in kwargs:
          self.parent       = kwargs["experiment"] # parent must be an experiment object
          self.totalDim     = self.parent.distribution.getDimension()
          self.xWeight = np.zeros(self.totalDim)
        elif "totalDim" in kwargs:
          self.parent     = None
          self.totalDim   = kwargs["totalDim"]
          self.xWeight = np.zeros(self.totalDim)
        else:
          ErrorMsg('The arguments must contain ("ID", "experiment" and "verbosity") or ("ID", "totalDim" and "verbosity").', fileName)
        self.verbosity    = kwargs["verbosity"]
        willExit          = False
      except:
        willExit = True
    if willExit:
      ErrorMsg('The arguments must contain ("block") or ("block" and "experiment") or ("ID", "totalDim" and "verbosity") or ("ID", "experiment" and "verbosity") or ("empty=True" and "verbosity").', fileName)


  def Build(self, xSample, xEval):
    # xSample must be of shape (nSample, nDimension)
    if xSample.shape[1] != self.totalDim:
      ErrorMsg('xSample[i,:] must be of length self.totalDim', fileName)

    if len(xSample) == len(xEval):
      prt('', 'green', self.verbosity)
      prt('Computing active direction ' + self.ID, 'green', self.verbosity)
      model = LinearRegression().fit(xSample, xEval)
      self.determinationCoeff = model.score(xSample, xEval)
      weights = model.coef_
      norm = norm2(weights)
      self.xWeight[:] = weights[:]/norm
      prt('Active direction ' + self.ID + ' computed', 'green', self.verbosity)
      prt('Determination coefficient = ' + str(self.determinationCoeff), 'green', self.verbosity)
    else:
      ErrorMsg('xSample[:,i] and xEval must have the same length', fileName)


  def _InitWithBlock(self, BLOCK):
    prt('', 'green', self.verbosity)
    prt('Initializing the Active direction from the input file.', 'green', self.verbosity)
    keys = BLOCK.GetKeysForBlock()
    if 'LoadFile' in keys:
      myFile = BLOCK.GetCharacterFromKey('LoadFile')
      if 'LoadPath' in keys:
        myPath = BLOCK.GetCharacterFromKey('LoadPath')
      else:
        myPath = ''
      self.Load(myFile, prefix=myPath)
    else:
      self.ID         = BLOCK.GetCharacterFromKey('ID')
      self.verbosity  = BLOCK.GetCharacterFromKey('outputLevel')


  def ComputeReducedFromTotal(self, totalVar):
    if totalVar.ndim == 1:
      redVar = 0
      if self.totalDim == len(totalVar):
        for j in range(self.totalDim):
          redVar += self.xWeight[j] * totalVar[j]
        return redVar
      else:
        ErrorMsg('The dimension of totalVar is not equal to the total dimension of the reduction method', fileName)
    elif totalVar.ndim == 2:
      nSample = totalVar.shape[0]
      sampleDim = totalVar.shape[1]
      if self.totalDim == sampleDim:
        xRedVar = np.zeros(nSample)
        for i in range(nSample):
          for j in range(self.totalDim):
            xRedVar[i] += self.xWeight[j] * totalVar[i,j]
        return xRedVar
      else:
        ErrorMsg('The dimension of totalVar[i] is not equal to the total dimension of the reduction method', fileName)
    else:
      ErrorMsg('Incompatible shape of totalVar', fileName)


  def ComputeTotalFromReduced(self, reducedVar):
    ErrorMsg('I do not know what TODO here because the problem is under-determined', fileName)


  def Save(self):
    self.savedItems['ID']                 = self.ID
    self.savedItems['TotalDim']           = self.totalDim
    self.savedItems['ReducedDim']         = self.reducedDim
    self.savedItems['ObjectType']         = self.objectType
    self.savedItems['xWeight']            = self.xWeight
    self.savedItems['DeterminationCoeff'] = self.determinationCoeff


  def Load(self, myInput, prefix=''):
    if isinstance(myInput, dict):
      self._LoadFromDictionary(myInput)
    elif isinstance(myInput, str):
      self._LoadFromFile(myInput, prefix=prefix)
    else:
      ErrorMsg('The argument "myInput" must be of type "dict" or "str" (name of file)', fileName)


  def _LoadFromFile(self, myFile, prefix=''):
    prt('', 'green', self.verbosity)
    prt('Loading the Active direction from the input file ' + myFile +'.', 'green', self.verbosity)
    myDict = DictFromH5(myFile, prefix=prefix)
    self._LoadFromDictionary(myDict)


  def _LoadFromDictionary(self, myDict):
    if myDict['ObjectType'] != self.objectType:
      ErrorMsg('The provided dictionnary defines an object of type ' + str(myDict['ObjectType']) + ' while it should define an object of type ' + self.objectType, fileName)
    self.ID                 = myDict['ID']
    self.totalDim           = myDict['TotalDim']
    self.reducedDim         = myDict['ReducedDim']
    self.xWeight            = myDict['xWeight']
    self.determinationCoeff = myDict['DeterminationCoeff']

    prt('ActiveDirection ' + self.ID + ' successfully loaded.', 'green', self.verbosity)
