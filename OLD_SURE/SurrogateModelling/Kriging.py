import h5py
import numpy as np
import sys
from Utils.I_O import prt
from Utils.I_O import ErrorMsg
from Utils.I_O import DictFromH5
import openturns as ot
import pickle
import random
import math
from SurrogateModelling.SurrogateModel import SurrogateModel as SM
from Numerics.Optimisation.Dichotomy  import DichotomyOptim  as DOpti
import sys

fileName = 'SurrogateModelling/Kriging.py'

class Kriging(SM):
  
  ######################
  #     Properties     #
  ######################
  multivariateBasis         = None
  covarianceModel           = None
  noiseLevel                = None
  result                    = None
  responseSurface           = None
  distribution              = None
  dimension                 = 0
  covarianceOptimization    = False
  noiseOptimization         = False

  trainingRatio             = 0
  trainingLength            = 0
  trainingEvals             = None
  trainingSamples           = None
  testSamples               = None
  testEvals                 = None
  algo                      = None

  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    willExit  = True
    self.savedItems   = {}
    self.objectType   = 'Surrogate_Kriging'
    self.noiseOptimization = False
    self.covarianceOptimization = False
    if "block" in kwargs and "experiment" in kwargs: # kwargs["block"] must be a block of "InputHandler"
      self.parent       = kwargs["experiment"] # Must be a valid experiment object with valid variables and distribution
      self._InitWithBlock(kwargs["block"])
      willExit          = False
    elif "empty" in kwargs:
      self.verbosity    = kwargs["verbosity"]
      willExit          = False
    else:
      try:
        self.ID           = kwargs["ID"]
        self.verbosity    = kwargs["verbosity"]
        if "experiment" in kwargs:
          self.parent       = kwargs["experiment"]
          self.distribution = self.parent.distribution
        elif "distribution" in kwargs:
          self.distribution = kwargs["distribution"]
        self.xVariable    = []          # List of variable objects TODO: Remove this ?
        willExit          = False
      except:
        willExit = True
      try:
        self.dimension    = self.distribution.getDimension()
      except:
        pass
    if willExit:
      ErrorMsg('The arguments must contain ("block" and "experiment") or ("ID", "experiment" and "verbosity") or ("ID", "distribution" and "verbosity") or ("empty=True" and "verbosity").', fileName)


  def Load(self, myInput, prefix=''):
    if isinstance(myInput, dict):
      self._LoadFromDictionary(myInput)
    elif isinstance(myInput, str):
      self._LoadFromFile(myInput, prefix=prefix)
    else:
      ErrorMsg('The argument "myInput" must be of type "dict" or "str" (name of file)', fileName)


  def _LoadFromFile(self, myFile, prefix=''):
    prt('', 'green', self.verbosity)
    prt('Loading the Kriging from the input file ' + myFile +'.', 'green', self.verbosity)
    myDict = DictFromH5(myFile, prefix=prefix)
    self._LoadFromDictionary(myDict)


  def _LoadFromDictionary(self, myDict):
    if myDict['ObjectType'] != self.objectType:
      ErrorMsg('The provided dictionary defines an object of type ' + str(myDict['ObjectType']) + ' while it should define an object of type ' + self.objectType, fileName)
    self.ID = myDict['ID']
    if myDict['LightSave']:
      self.result           = None
      self.responseSurface  = pickle.loads(myDict['ResponseSurface'])
    else:
      self.result           = pickle.loads(myDict['Result'])
      self.responseSurface  = self.result.getMetaModel()
      self.dimension        = myDict['Dimension']
    prt('Kriging ' + self.ID + ' successfully loaded.', 'green', self.verbosity)


  def Copy(self, other):
    pass


  def Save(self, light = False):
    self.savedItems['ID']                = self.ID
    self.savedItems['ObjectType']        = self.objectType
    self.savedItems['LightSave']         = light
    self.savedItems['Dimension']         = self.dimension
    if light:
      self.savedItems['ResponseSurface'] = pickle.dumps(self.responseSurface)
    else:
      self.savedItems['Result']          = pickle.dumps(self.result)


  def Evaluate(self, inputs):
    return self.responseSurface(inputs) #TODO: accept arrays of inputs like the 1D Spline


  ################################################################
  # Building the optimal noisy kriging surrogate
  ############################################################
  #-----------------------------------------------------------------
  # Defining the error function:
  #-----------------------------------------------------------------
  def _ErrorEvaluation(self, rSurf):
    nTest         = len(self.testEvals)
    modelEvals    = np.zeros(nTest)
    for i in range(nTest):
      modelEvals[i] = rSurf(self.testSamples[i, 0])
    error         = 0.0
    for i in range(nTest):
      error += (modelEvals[i]-self.testEvals[i])*(modelEvals[i]-self.testEvals[i])
    #error /= (nTest*abs(np.mean(modelEvals)))
    return error
  #-----------------------------------------------------------------
  # Function to build surrogate and evaluate error:
  #-----------------------------------------------------------------
  def _objectiveFunc(self, noiseLevel):
    self.DefineNoise(noiseLevel)
    krig, rSurf = self.Build_single()
    error       = self._ErrorEvaluation(rSurf)
    return error
  #-----------------------------------------------------------------
  # Surrogate building:
  #-----------------------------------------------------------------
  def Build(self, xSample, xEval, isBootstrap=False, bootMultiplier=10,
            trainingRatio=0.8, fixSeed=False, lowerBound=1e-10,
            upperBound=1e0, nInitialSamples=10, relativeTolerance=1.0e-3):

    if len(xSample) != len(xEval):
      ErrorMsg('Samples and evaluations must have the same size.', fileName)
    self.trainingRatio = trainingRatio

    if not self.noiseOptimization:                # Noise level is set
      self.trainingSamples  = xSample
      self.trainingEvals    = xEval
      prt('', 'green', self.verbosity)
      prt('Building Kriging ' + self.ID + '...', 'green', self.verbosity)
      self.result, self.responseSurface = self.Build_single()
      prt('Kriging ' + self.ID + ' computed', 'green', self.verbosity)

    else:                                         # Noise level must be optimized
      ErrorMsg('Noise level optimization is not functional yet. Please set the noise level manually.', fileName)
      # Step 1: Find the optimal noise level # TODO: Keep the best surrogate in memory to avoid rebuilding it
      prt("", 'green', self.verbosity)
      prt('Dividing sample set...', 'green', self.verbosity)
      self._DivideSamples(xSample, xEval, isBootstrap, bootMultiplier, None, fixSeed)
      prt('Finding optimal input noise level...', 'green', self.verbosity)
      optimalNoiseLevel = DOpti(self._objectiveFunc, lowerBound, upperBound, nInitialSamples, relativeTolerance,
                                dichoType='geometric')

      # Step 2: Build the optimal Kriging
      dim = len(self.trainingSamples[0,:])
      dataSize = len(self.trainingEvals)
      samples = ot.Sample(dataSize, dim)
      samples[:,:] = self.trainingSamples[:,:]
      evals = ot.Sample(dataSize, 1)
      for i in range(dataSize):
        evals[i, 0] = self.trainingEvals[i]

      algo = ot.KrigingAlgorithm(samples, evals, self.covarianceModel, self.multivariateBasis, True)
      algo.setNoise([optimalNoiseLevel] * len(xEval))
      startingPoint = ot.LHSExperiment(ot.ComposedDistribution([ot.LogUniform(1e-1, 1e2)]*self.dimension), 50).generate()
      algo.setOptimizationAlgorithm(ot.MultiStart(ot.TNC(), startingPoint))
      algo.setOptimizationBounds(ot.Interval([0.1]*self.dimension, [1e2]*self.dimension))
      algo.setOptimizeParameters(self.covarianceOptimization)
      algo.run()
      self.result           = algo.getResult()
      self.responseSurface  = self.result.getMetaModel()
      prt("Kriging " + self.ID + " computed. Noise level = " + str(optimalNoiseLevel), 'green', self.verbosity)

  #-----------------------------------------------------------------
  # Single Kriging computation:
  #-----------------------------------------------------------------
  def Build_single(self):
    xSample = self.trainingSamples
    xEval   = self.trainingEvals
    dim = len(xSample[0,:])
    dataSize = len(xEval)
    samples = ot.Sample(dataSize, dim)
    samples[:,:] = xSample[:,:]
    evals = ot.Sample(dataSize, 1)
    for i in range(dataSize):
      evals[i, 0] = xEval[i]

    algo = ot.KrigingAlgorithm(samples, evals, self.covarianceModel, self.multivariateBasis, True)
    if self.noiseLevel is not None:
      algo.setNoise([self.noiseLevel] * len(xEval))
    startingPoint = ot.LHSExperiment(ot.ComposedDistribution([ot.LogUniform(1e-2, 1e4)] * self.dimension), 500).generate()
    algo.setOptimizationAlgorithm(ot.MultiStart(ot.TNC(), startingPoint))
    algo.setOptimizationBounds(ot.Interval([1e-2]*self.dimension, [1e4]*self.dimension))
    algo.setOptimizeParameters(self.covarianceOptimization)
    algo.run()
    result           = algo.getResult()
    responseSurface  = result.getMetaModel()
    return result, responseSurface


  ################################################################
  # Managing the sample sets
  ############################################################
  def _ClearSampleSets(self):
    self.trainingEvals = None
    self.trainingSamples = None
    self.testSamples = None
    self.testEvals = None

  def _DivideSamples(self, samples, evaluations, isBootstrap, bootMultiplier, weights, fixSeed):
    # Initialisation
    self._ClearSampleSets()
    nSample = len(samples)

    if fixSeed:
      random.seed(0)
    if not isBootstrap:
      self.trainingLength = int(math.floor(self.trainingRatio * nSample))
      testIndexList = np.array(range(nSample))
      trainingIndexList = np.zeros(self.trainingLength, dtype=int)
      for i in range(len(trainingIndexList)):
        index = random.sample(testIndexList, 1)
        trainingIndexList[i] = index[0]
        testIndexList = testIndexList[testIndexList != index[0]]

      # Training set:
      self.trainingSamples = samples[trainingIndexList]
      self.trainingEvals = evaluations[trainingIndexList]
      # Test set:
      self.testSamples = samples[testIndexList]
      self.testEvals = evaluations[testIndexList]
    else:
      ErrorMsg('The bootstrap sampling method is not yet implemented. Please try again later.', fileName)








  ################################################################
  # Defining Kriging options
  ############################################################
  def DefineTrend(self, trend):
    if trend == 'Constant':
      self.multivariateBasis = ot.ConstantBasisFactory(self.dimension).build()
    elif trend == 'Linear':
      self.multivariateBasis = ot.LinearBasisFactory(self.dimension).build()
    elif trend == 'Quadratic':
      self.multivariateBasis = ot.QuadraticBasisFactory(self.dimension).build()
    elif trend == 'PCE':
      #self.multivariateBasis = ot.ConstantBasisFactory(self.dimension).build()
      ErrorMsg('Trend type "PCE" has not been implemented yet.', fileName)
    else:
      ErrorMsg('Trend type "' + trend +'" is not valid. Valid types are "Constant", "Linear", Quadratic" and "PCE"', fileName)

  def DefineCovariance(self, covarianceType):
    if covarianceType == 'Matern':
      self.covarianceModel = ot.MaternModel(self.dimension)
    else:
      ErrorMsg('Covariance type "' + covarianceType +'" is not supported for now.', fileName)

  def SetCovarianceOptimization(self, bool=True):
    self.covarianceOptimization = bool

  def SetNoiseOptimization(self, bool=True):
    self.noiseOptimization = bool

  def DefineNoise(self, noiseLevel):
    self.noiseLevel = noiseLevel

  def SetDistribution(self, distribution):
    self.distribution = distribution

  ################################################################
  # Block initialisation
  ############################################################
  def _InitWithBlock(self, BLOCK):
    prt('', 'green', self.verbosity)
    prt('Initializing the Kriging from the input file.', 'green', self.verbosity)
    
    self.distribution = self.parent.distribution
    self.dimension = self.distribution.getDimension()
    # ----- General information -----
    keys = BLOCK.GetKeysForBlock()
    self.ID         = BLOCK.GetCharacterFromKey('ID')
    self.verbosity  = BLOCK.GetCharacterFromKey('outputLevel')
    trendType       = BLOCK.GetCharacterFromKey('trendType')
    covarianceModel = BLOCK.GetCharacterFromKey('covarianceModel')
    self.DefineTrend(trendType)
    self.DefineCovariance(covarianceModel)
    self.SetCovarianceOptimization()
    if "optimizeNoiseLevel" in keys:
      self.SetNoiseOptimization(BLOCK.GetLogicalFromKey('optimizeNoiseLevel'))
    elif "noiseLevel" in keys:
      self.DefineNoise(BLOCK.GetRealFromKey('noiseLevel'))


