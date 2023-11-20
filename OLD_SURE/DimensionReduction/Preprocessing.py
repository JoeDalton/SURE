import numpy as np
import math
import sys
from Utils.I_O import prt, DumpObject, ErrorMsg
from copy import copy as cp

fileName = 'DimensionReduction/Preprocessing.py'

class Preprocessing:
  ######################
  #     Properties     #
  ######################
  parent             = None
  savedItems         = {}
  verbosity          = -1
  preType            = ''
  ID                 = ''
  isGlobal           = False

  ###################
  #     Methods     #
  ###################
  def __init__(self, ID, parent=None, verbosity=-1):
    self.parent     = parent
    self.verbosity  = verbosity
    self.savedItems = {}
    self.totalDim   = 0
    self.reducedDim = 0
    self.preType    = ''
    self.ID         = ID
    self.isGlobal   = False

  def Dump(self, directory='./'):
    DumpObject(self.savedItems, directory)
    prt(self.redType + "_" + self.ID + " : Dumped", "green", self.verbosity)



class CenterScaling(Preprocessing):
  # First centers, then scales
  ######################
  #     Properties     #
  ######################
  myCentering = None
  myScaling   = None
  
  ###################
  #     Methods     #
  ###################
  def __init__(self, ID, data, scaleType, isCGlobal=False, isSGlobal = False, parent=None, verbosity=-1):
    # "data" can be a numpy vector or a numpy array (nFeature, nSample)
    # Initialize
    self.ID                = ID 
    self.parent            = parent
    self.verbosity         = verbosity
    self.savedItems        = {}
    self.preType           = 'CenterScaling'
    self.myScaling         = Scaling(ID + '_Scaling', data, scaleType, isGlobal=isSGlobal, parent=self, verbosity=self.verbosity-1)
    print("scaling ok")    
    self.myCentering       = Centering(ID + '_Centering', data, isGlobal=isCGlobal, parent=self, verbosity=self.verbosity-1)
    print("centering ok")    
    self.isGlobal          = isCGlobal and isSGlobal
     
  def Save(self):
    self.savedItems['ID']                 = self.ID
    self.savedItems['Verbosity']          = self.verbosity
    self.savedItems['ObjectType']         = "Preprocessing_" + self.preType
    self.savedItems['xMean']              = self.xMean

  def ComputePro(self, orig):
    processed = self.myScaling.ComputePro(self.myCentering.ComputePro(orig))
    return processed

  def ComputeOrig(self, pro):
    orig = self.myCentering.ComputeOrig(self.myScaling.ComputeOrig(pro))
    return orig


class Centering(Preprocessing):
  ######################
  #     Properties     #
  ######################
  xMean            = []
  
  ###################
  #     Methods     #
  ###################
  def __init__(self, ID, data, isGlobal=False, parent=None, verbosity=-1):
    # "data" can be a numpy vector or a numpy array (nFeature, nSample)
    # Initialize
    self.ID                = ID 
    self.parent            = parent
    self.verbosity         = verbosity
    self.savedItems        = {}
    self.preType           = 'Centering'
    self.isGlobal          = isGlobal
    # Compute the means of the data
    self.xMean = []
    if isinstance(data, np.ndarray):
      if len(data.shape) == 1:
        self.xMean.append(np.mean(data))
      else:
        for i in range(data.shape[0]):
          if self.isGlobal:
            self.xMean.append(np.mean(data))
          else:
            self.xMean.append(np.mean(data[i,:]))
    else:
      ErrorMsg('The "data" input must be a np.array', fileName)

  def Save(self):
    self.savedItems['ID']                 = self.ID
    self.savedItems['Verbosity']          = self.verbosity
    self.savedItems['ObjectType']         = "Preprocessing_" + self.preType
    self.savedItems['xMean']              = self.xMean

  def ComputePro(self, orig):
    pro = np.zeros_like(orig)
    if isinstance(orig, np.ndarray):
      if len(orig.shape) == 1:
        pro[:] = orig[:] - self.xMean[:]
      else:
        if orig.shape[0] == len(self.xMean):
          for i in range(orig.shape[0]):
            pro[i,:] = orig[i,:] - self.xMean[i]
        else:
          ErrorMsg('The "orig" input must be of shape (nFeature,nSample)', fileName)
    else:
      ErrorMsg('The "orig" input must be a numpy array of size (nSample,) (single feature) or (nFeature,nSample)', fileName)
    return pro
  
  def ComputeOrig(self, pro):
    orig = np.zeros_like(pro)
    if isinstance(pro, np.ndarray):
      if len(pro.shape) == 1:
        orig[:] = pro[:] + self.xMean[:]
      else:
        if pro.shape[0] == len(self.xMean):
          for i in range(pro.shape[0]):
            orig[i,:] = pro[i,:] + self.xMean[i]
        else:
          ErrorMsg('The "pro" input must be of shape (nFeature,nSample)', fileName)
    else:
      ErrorMsg('The "pro" input must be a numpy array of size (nSample,) (single feature) or (nFeature,nSample)', fileName)
    return orig


class Scaling(Preprocessing):
  ######################
  #     Properties     #
  ######################
  xScale            = []
  method            = ''
  
  ###################
  #     Methods     #
  ###################
  def __init__(self, ID, data, scaleType, isGlobal=False, parent=None, verbosity=-1):
    # "data" can be a numpy vector or a numpy array (nFeature, nSample) #TODO: transpose the formulation of the preprocessing ? This is not ideal and may cause bugs later on
    # Initialize
    self.ID                = ID 
    self.parent            = parent
    self.verbosity         = verbosity
    self.savedItems        = {}
    self.preType           = 'Scaling_' + scaleType
    self.method            = scaleType
    self.isGlobal          = isGlobal
    # Compute the means of the data
    self.xScale = []
    if isinstance(data, np.ndarray):
      if len(data.shape) == 1:
        self.xScale.append(self._ComputeVectorScale(data, scaleType))
      else:
        for i in range(data.shape[0]):
          if self.isGlobal:
            self.xScale.append(self._ComputeVectorScale(data, scaleType))
          else:
            self.xScale.append(self._ComputeVectorScale(data[i,:], scaleType))
    else:
      ErrorMsg('The "data" input must be a np.array', fileName)

  def _ComputeVectorScale(self, vector, scaleType):
    #print(vector)
    #print(len(vector))
    if scaleType == 'AUTO':
      # All samples become equally important: this is the "standard score"
      scale = np.std(vector)
    elif scaleType == 'RANGE':
      # All samples become equally important
      scale = (vector.max()-vector.min())
    elif scaleType == 'PARETO':
      # Reduces the importance of large values, but keep data structure partially intact
      scale = math.sqrt(np.std(vector))
    elif scaleType == 'VAST':
      # Focuses on small fluctuations, aims at robustness
      scale = np.mean(vector) / (np.std(vector)**2)
    elif scaleType == 'LEVEL':
      # Focuses on relative response
      scale = np.mean(vector)
    elif scaleType == 'MAX':
      scale = np.max(vector)
    else:
      ErrorMsg('The "scaleType" input must be "AUTO", "RANGE", "PARETO", "VAST", "LEVEL" or "MAX"', fileName)
    if scale == 0:
      ErrorMsg('The scaling type chosen is incompatible with your data. (scale = 0). Please choose another type of scaling.', fileName)
    return scale

  def Save(self):
    self.savedItems['ID']                 = self.ID
    self.savedItems['Verbosity']          = self.verbosity
    self.savedItems['ObjectType']         = "Preprocessing_" + self.preType
    self.savedItems['xScale']             = self.xScale
    
  def ComputePro(self, orig):
    pro = np.zeros_like(orig)
    if isinstance(orig, np.ndarray):
      if len(orig.shape) == 1:
        pro[:] = orig[:] / self.xScale[:]
      else:
        if orig.shape[0] == len(self.xScale):
          for i in range(orig.shape[0]):
            pro[i,:] = orig[i,:] / self.xScale[i]
        else:
          ErrorMsg('The "orig" input must be of shape (nFeature,nSample)', fileName)
    else:
      ErrorMsg('The "orig" input must be a numpy array of size (nSample,) (single feature) or (nFeature,nSample)', fileName)
    return pro
  
  def ComputeOrig(self, pro):
    orig = np.zeros_like(pro)
    if isinstance(pro, np.ndarray):
      if len(pro.shape) == 1:
        orig[:] = pro[:] * self.xScale[:]
      else:
        if pro.shape[0] == len(self.xScale):
          for i in range(pro.shape[0]):
            orig[i,:] = pro[i,:] * self.xScale[i]
        else:
          ErrorMsg('The "pro" input must be of shape (nFeature,nSample)', fileName)
    else:
      ErrorMsg('The "pro" input must be a numpy array of size (nSample,) (single feature) or (nFeature,nSample)', fileName)
    return orig


########################################################
"""
This is a test for a generic preprocessor
TODO: Test this extensively, because I fear it has be done a bit too quickly
TODO: Move this inside a "Preprocessor" object ?
"""

def Preprocessor(xSample, xEval, BLOCK):
  from ActiveDirection import ActiveDirection
  # Determine preprocessings to do
  keys        = BLOCK.GetKeysForBlock()
  # Pre-processing of sampling points
  if "xInPreProc" in keys:
    temp = BLOCK.GetCharacterFromKey('xInPreProc')
    if isinstance(temp, list):
      xInPreProc  = temp
    else:
      xInPreProc  = [temp]
  else:
    xInPreProc    = []
  # Pre-Processing of sampled QoIs
  if "xOutPreProc" in keys:
    temp = BLOCK.GetCharacterFromKey('xOutPreProc')
    if isinstance(temp, list):
      xOutPreProc = temp
    else:
      xOutPreProc = [temp]
  else:
    xOutPreProc   = []

  # Preprocess the QoIs # Must be done first so that an eventual active subspace will be correctly computed
  if xEval is None:
    xProcEval = None
  else:
    tempOut = cp(xEval)
    for proc in xOutPreProc:
      BLOCK.GetSubBlock(proc)
      preType       = BLOCK.GetCharacterFromKey('Type')
      if preType == 'LogTransform':
        tempOut = np.log(tempOut)
      else:
        ErrorMsg('The supported QoI preprocessing are: "LogTransform".', fileName) # TODO: Accept more, for example the scalings, and perhaps also other analyticals if necessary
    xProcEval = tempOut

  # Preprocess the sample points
  if xSample is None:
    xProcSample = None
  else:
    tempIn = cp(xSample)
    for proc in xInPreProc:
      BLOCK.GetSubBlock(proc)
      preType       = BLOCK.GetCharacterFromKey('Type')
      if preType == 'ActiveDirection':
        #saveFile      = BLOCK.GetCharacterFromKey('SaveFile')
        #savePath      = BLOCK.GetCharacterFromKey('SavePath')
        myAD = ActiveDirection(block=BLOCK) # Warning: not rigorous if it is not the first applied...
        #myAD.Build(xSample, xEval) # TODO: Add the option somewhere to build and save the AD ?
        #myAD.Save()
        #myAD.Dump(fileName=saveFile, path=savePath)
        tempIn = myAD.ComputeReducedFromTotal(tempIn)
        # DEBUG: this trick is necessary to work with kriging. TODO: Try to improve ?
        newTempIn = np.zeros((len(tempIn), 1))
        for i in range(len(tempIn)):
          newTempIn[i, 0] = tempIn[i]
        tempIn = newTempIn
      elif preType == 'LogTransform':
        tempIn = np.log(tempIn)
      else:
        ErrorMsg('The supported sample points preprocessing are "ActiveDirection" and "LogTransform".', fileName) # TODO: Accept more, for example the scalings, and perhaps also other analyticals if necessary
    xProcSample = tempIn

  return xProcSample, xProcEval
