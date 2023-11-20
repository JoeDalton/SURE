import sys
import math
import random
import pickle
import numpy                                                    as np
from Utils.I_O                        import prt
from Utils.I_O                        import ErrorMsg
from Utils.I_O                        import DictFromH5
from SurrogateModel                   import SurrogateModel
from scipy.interpolate                import UnivariateSpline
from sklearn.linear_model             import LinearRegression
from Numerics.Optimisation.Dichotomy  import DichotomyOptim     as DOpti

fileName = 'SurrogateModelling/OneDSpline.py'

class OneDSpline(SurrogateModel):
  
  ######################
  #     Properties     #
  ######################
  trainingRatio             = 0.0
  responseSurface           = None
  smoothingFactor           = 0.0
  nKnots                    = 0
  trainingLength            = 0
  trainingEvals             = None 
  trainingWeights           = None 
  trainingSamples           = None 
  testSamples               = None 
  testEvals                 = None 
  testWeights               = None
  dimension                 = 1
  distribution              = None
  


  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    self.savedItems   = {}
    self.dimension    = 1
    self.distribution = None
    willExit = True
    self.objectType   = 'Surrogate_1DSpline'
    if "block" in kwargs and "experiment" in kwargs: # kwargs["block"] must be a block of "InputHandler"
      self.parent       = kwargs["experiment"] # Must be a valid experiment object with valid variables and distribution
      self._InitWithBlock(kwargs["block"]) 
      willExit          = False
    elif "empty" in kwargs and "verbosity" in kwargs:
      self.verbosity    = kwargs["verbosity"]
      willExit          = False
    else:
      try:
        self.ID           = kwargs["ID"]
        self.verbosity    = kwargs["verbosity"]
        if "experiment" in kwargs:
          self.parent       = kwargs["experiment"]
        willExit          = False
      except:
        willExit = True
    if willExit:
      ErrorMsg('The arguments must contain ("block" and "experiment") or ("ID", "experiment" and "verbosity") or ("ID" and "verbosity") or ("empty=True" and "verbosity").', fileName)


  def Load(self, myInput, prefix=''):
    if isinstance(myInput, dict):
      self._LoadFromDictionary(myInput)
    elif isinstance(myInput, str):
      self._LoadFromFile(myInput, prefix=prefix)
    else:
      ErrorMsg('The argument "myInput" must be of type "dict" or "str" (name of file)', fileName)


  def _LoadFromFile(self, myFile, prefix=''):
    prt('', 'green', self.verbosity)
    prt('Loading the OneDSpline from the input file ' + myFile +'.', 'green', self.verbosity)
    myDict = DictFromH5(myFile, prefix=prefix)
    self._LoadFromDictionary(myDict)


  def _LoadFromDictionary(self, myDict):
    if myDict['ObjectType'] != self.objectType:
      ErrorMsg('The provided dictionary defines an object of type ' + str(myDict['ObjectType']) + ' while it should define an object of type ' + self.objectType, fileName)
    self.ID = myDict['ID']
    self.responseSurface  = pickle.loads(myDict['ResponseSurface'])
    prt('OneDSpline ' + self.ID + ' successfully loaded.', 'green', self.verbosity)


  def Copy(self, other):
    pass


  def Save(self):
    self.savedItems['ID']              = self.ID
    self.savedItems['ObjectType']      = self.objectType
    self.savedItems['ResponseSurface'] = pickle.dumps(self.responseSurface)


  def Evaluate(self, inputs): # inputs can be either an input scalar or a vector of vector one which the spline must be evaluated
    if isinstance(inputs, list) or isinstance(inputs, np.ndarray):
      return self.responseSurface(inputs)
    else:
      return float(self.responseSurface(inputs))

  def SetDistribution(self, distribution):
    self.distribution = distribution

  ################################################################
  # Managing the sample sets
  ############################################################
  def _DivideSamples(self, samples, evaluations, isBootstrap, bootMultiplier, weights, fixSeed):
    # Initialisation
    self._ClearSampleSets()
    nSample = len(samples)
    
    #if weights is None: # DEBUG: Do we really wish to see non-unity weights for this ?
    weights = np.ones(nSample)
    if fixSeed:
      random.seed(0)
    if not isBootstrap:
      self.trainingLength = int(math.floor(self.trainingRatio * nSample))
      testIndexList       = np.array(range(nSample))
      trainingIndexList   = np.zeros(self.trainingLength, dtype=int)
      for i in range(len(trainingIndexList)):
        index = random.sample(testIndexList, 1)
        trainingIndexList[i] = index[0]
        testIndexList = testIndexList[testIndexList != index[0]]
 
      # Training set:
      trainingSamples = samples[trainingIndexList]
      trainingEvals   = evaluations[trainingIndexList]
      trainingWeights = weights[trainingIndexList]
      # Sort training set for spline fitting:
      trainSort       = trainingSamples.argsort()
      self.trainingEvals   = trainingEvals[trainSort]
      self.trainingWeights = trainingWeights[trainSort]
      self.trainingSamples = trainingSamples[trainSort]
      # Test set:
      self.testSamples     = samples[testIndexList]
      self.testEvals       = evaluations[testIndexList]
      self.testWeights     = weights[testIndexList]
    else:
      ErrorMsg('The bootstrap sampling method is not yet implemented. Please try again later.', fileName)


  def _ClearSampleSets(self):
    self.trainingEvals   = None 
    self.trainingWeights = None 
    self.trainingSamples = None 
    self.testSamples     = None 
    self.testEvals       = None 
    self.testWeights     = None 



  ################################################################
  # Building the optimal spline surrogate
  ############################################################
  #-----------------------------------------------------------------
  # Defining the error function:
  #-----------------------------------------------------------------
  def _ErrorEvaluation(self, spl):
    maxNKnots     = 10
    nTest         = len(self.testSamples)
    modelEvals    = spl(self.testSamples)
    error         = 0.0
    for i in range(nTest):
      error += (modelEvals[i]-self.testEvals[i])*(modelEvals[i]-self.testEvals[i])*self.testWeights[i]
    error /= (nTest*abs(np.mean(modelEvals)))
    nKnots = len(spl.get_knots())
    if nKnots > maxNKnots:
      error += (nKnots-maxNKnots)
    return error
  #-----------------------------------------------------------------
  # Function to compute spline and evaluate error:
  #-----------------------------------------------------------------
  def _objectiveFunc(self, smoothingFactor):
    smoothing = smoothingFactor * self.trainingLength
    spl       = UnivariateSpline(self.trainingSamples, self.trainingEvals, w=self.trainingWeights, k=3, s=smoothing)
    error     = self._ErrorEvaluation(spl)
    return error
  #-----------------------------------------------------------------
  # Launching the optimisation:
  #-----------------------------------------------------------------
  def Build(self, samples=None, evaluations=None, isBootstrap=False, bootMultiplier=10, 
            trainingRatio=0.8, weights=None, fixSeed=False, lowerBound=0.001, 
            upperBound=10, nInitialSamples=10, relativeTolerance=1.0e-3):
    

    # Initialisation
    self.trainingRatio  = trainingRatio
    #if samples is None and evaluations is None:
    #  samples     = self.parent.samples
    #  evaluations = self.parent.evaluations
    #if samples is None:
    #  prt('SurrogateModelling/1DSpline.py/Build: Please provide samples.', 'red', self.verbose)
    #if evaluations is None:
    #  prt('SurrogateModelling/1DSpline.py/Build: Please provide evaluations.', 'red', self.verbose)
    if len(samples) != len(evaluations):
      ErrorMsg('Samples and evaluations must have the same size.', fileName)

    # Building the training set and the test set      
    prt("", 'green', self.verbosity)
    prt("Computing spline surrogate " + self.ID, 'green', self.verbosity)
    self._DivideSamples(samples, evaluations, isBootstrap, bootMultiplier, weights, fixSeed)
   
    optimalSmoothingFactor = DOpti(self._objectiveFunc, lowerBound, upperBound, nInitialSamples, relativeTolerance, dichoType='geometric')

    absSmoothing           = optimalSmoothingFactor * self.trainingLength
    self.responseSurface   = UnivariateSpline(self.trainingSamples, self.trainingEvals, w=self.trainingWeights, k=3, s=absSmoothing)
    self.nKnots            = len(self.responseSurface.get_knots())
    prt("Spline surrogate " + self.ID + " computed. Number of knots = " + str(self.nKnots), 'green', self.verbosity)
    # Purging the training set and the test set
    self._ClearSampleSets()


  def _InitWithBlock(self, BLOCK):
    prt('', 'green', self.verbosity)
    prt('Initializing the OneDSpline from the input file.', 'green', self.verbosity)
    
    self.ID         = BLOCK.GetCharacterFromKey('ID')
    self.verbosity  = BLOCK.GetCharacterFromKey('outputLevel')
