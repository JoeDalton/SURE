import h5py
import numpy as np
import sys
from Utils.I_O import prt
from Utils.I_O import ErrorMsg, WarningMsg
from Utils.I_O import DictFromH5
import openturns as ot
import pickle
from SurrogateModel import SurrogateModel
import sys

fileName = 'SurrogateModelling/PCE.py'

class PCE(SurrogateModel):
  
  ######################
  #     Properties     #
  ######################
  multivariateBasis         = None
  hyperbolicQuasiNorm       = 1.0
  truncatureBasisStrategy   = None
  coeffEvaluationStrategy   = None
  polyCoeffs                = None
  polyIndices               = None
  result                    = None
  responseSurface           = None
  dimension                 = 0
  distribution              = None
  checkBasis                = False


  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    willExit  = True
    self.savedItems   = {}
    self.objectType   = 'Surrogate_PCE'
    if "block" in kwargs and "experiment" in kwargs: # kwargs["block"] must be a block of "InputHandler"
      self.parent       = kwargs["experiment"] # Must be a valid experiment object with valid variables and distribution
      self._InitWithBlock(kwargs["block"]) 
      willExit          = False
    elif "empty" in kwargs:
      self.verbosity    = kwargs["verbosity"]
      willExit          = False
    else:
      try:
        self.verbosity = kwargs["verbosity"]
        self.ID           = kwargs["ID"]
        if "checkBasis" in kwargs:
          self.checkBasis = kwargs["checkBasis"]
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
    prt('Loading the PCE from the input file ' + myFile +'.', 'green', self.verbosity)
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
      self.polyCoeffs       = self.result.getCoefficients()
      self.polyIndices      = self.result.getIndices()
      self.responseSurface  = self.result.getMetaModel()
      self.dimension        = myDict['Dimension']
    prt('PCE ' + self.ID + ' successfully loaded.', 'green', self.verbosity)


  def Copy(self, other):
    pass


  def Save(self, light = False):
    self.savedItems['ID']                = self.ID
    self.savedItems['ObjectType']        = self.objectType
    self.savedItems['LightSave']         = light
    if light:
      self.savedItems['ResponseSurface'] = pickle.dumps(self.responseSurface)
    else:
      self.savedItems['Result']          = pickle.dumps(self.result)
  

  def Build(self, xSample, xEval, **kwargs):
    prt('', 'green', self.verbosity)
    prt('Converting input data...', 'green', self.verbosity)
    dim = len(xSample[0,:])
    dataSize = len(xEval)
    if 'xWeight' in kwargs:
      xWeight = kwargs['xWeight'] * dataSize
    else:
      xWeight = np.ones(dataSize)
    samples = ot.Sample(dataSize, dim)
    samples[:,:] = xSample[:,:]
    evals = ot.Sample(dataSize, 1)
    for i in range(dataSize):             # For some reason the vectorized expression was not accepted here...
      evals[i, 0] = xEval[i] * xWeight[i] # Multiply every evaluation by its weight (1 for Monte-Carlo, whatever needed for cubature)

    prt('', 'green', self.verbosity)
    prt('Building PCE ' + self.ID + '...', 'green', self.verbosity)

    algo = ot.FunctionalChaosAlgorithm(samples, evals, self.distribution, self.truncatureBasisStrategy, self.coeffEvaluationStrategy)
    algo.run()
    self.result           = algo.getResult()
    self.polyCoeffs       = self.result.getCoefficients()
    self.polyIndices      = self.result.getIndices()
    self.responseSurface  = self.result.getMetaModel()
    prt('PCE ' + self.ID + ' computed', 'green', self.verbosity)


  def Evaluate(self, inputs): # inputs must be a vector containing the parameters of the PCE
    return self.responseSurface(inputs)


  def DefinePolyBasis(self, quasiNorm=1.0, weights=None):
    if weights is None:
      weights = [1.0] * self.dimension
    # Create the univariate polynomial family collection
    # which regroups the polynomial families for each direction
    polyColl = ot.PolynomialFamilyCollection(self.dimension)
    for dimIdx in range(self.dimension):
      # TODO: If necessary, add other basis polynomials for other distributions
      if self.distribution.getMarginal(dimIdx) == ot.Normal(0.0,1.0):
        polyColl[dimIdx] = ot.HermiteFactory()
      elif self.distribution.getMarginal(dimIdx) == ot.Uniform(0.0,1.0) or self.distribution.getMarginal(dimIdx) == ot.Uniform(-1.0,0.0):
        polyColl[dimIdx] = ot.LegendreFactory()
      else:
        ErrorMsg('As yet unsupported distribution for PCE: ' + str(self.distribution.getMarginal(dimIdx)), fileName)

    # Create the enumeration function and save parameters
    self.hyperbolicQuasiNorm = quasiNorm
    self.savedItems['Dimension']           = self.dimension
    self.savedItems['HyperbolicQuasiNorm'] = quasiNorm
    self.savedItems['AnisotropyWeights']   = weights
    if self.hyperbolicQuasiNorm == 1.0:
      # Linear enumeration (Basic): The basis polynomials are taken in increasing degree order
      enumerateFunction = ot.LinearEnumerateFunction(self.dimension)
    else:
      if weights==[1.0] * self.dimension:
        # Hyperbolic enumeration: Based on the principle of the sparsity-of-effects:
        # The variability is assumed to come mainly from single contributions and low-order interactions.
        # Hence, this option allows to dismiss the high-order interaction polynomials whose contributions will likely be close to zero
        enumerateFunction = ot.HyperbolicAnisotropicEnumerateFunction(self.dimension, self.hyperbolicQuasiNorm)
      else:
        # weights must contain dim elements
        # If weights are specified, an anisotropic enumeration is used.
        # Same principle as above, but uncertain variables with high weights will get more high-order terms
        enumerateFunction = ot.HyperbolicAnisotropicEnumerateFunction(weights, self.hyperbolicQuasiNorm)

    # Create the multivariate orthonormal basis
    # which is the the cartesian product of the univariate basis
    self.multivariateBasis = ot.OrthogonalProductPolynomialFactory(polyColl, enumerateFunction)


  def DefineTruncatureStrat(self, strategy="Fixed", **kwargs):
    self.savedItems['TruncatureStrategy'] = strategy
    if 'maxNumPolys' in kwargs:
      maxNumPolys = kwargs['maxNumPolys']
      if 'maxTotDegree' in kwargs:
        WarningMsg('In the presence of the argument "maxNumPolys", "maxTotDegree" will be ignored', fileName, self.verbosity)
    elif 'maxTotDegree' in kwargs:
      maxNumPolys = self.multivariateBasis.getEnumerateFunction().getStrataCumulatedCardinal(kwargs['maxTotDegree'])
    else:
      ErrorMsg('Either argument "maxNumPolys" or "maxTotDegree" must be provided', fileName)
    self.savedItems['MaxNumPolys']        = maxNumPolys
    if strategy == "Fixed":
      # The first maxNumPolys polynomials of the basis
      self.truncatureBasisStrategy = ot.FixedStrategy(self.multivariateBasis, maxNumPolys)
    elif strategy == "Sequential":
      # Among the maxNumPolys first polynomials
      # of the multivariate basis, keep those verifying the convergence criterion
      self.truncatureBasisStrategy = ot.SequentialStrategy(self.multivariateBasis, maxNumPolys)
      WarningMsg('PCE.py/DefineTruncatureStrat: The robustness of the "Sequential" strategy has not yet been verified. Use at your own risk.', fileName, self.verbosity)
    elif strategy == "Cleaning":
      if "mostSignificant" in kwargs:
        mostSignificant = kwargs["mostSignificant"]
      else:
        print("Error")
        quit()
      if "significanceFactor" in kwargs:
        significanceFactor = kwargs["significanceFactor"]
      else:
        print("Error")
        quit()
      self.savedItems['MostSignificantPolys'] = mostSignificant
      self.savedItems['SignificantFactor']    = significanceFactor
      # Among the maxNumPolys first polynomials,
      # Keep those which have the mostSignificant most significant contributions
      # with significance criterion significanceFactor 
      # The last argument indicates if we are interested in the online monitoring of the current basis update
      # (removed or added coefficients)
      self.truncatureBasisStrategy = ot.CleaningStrategy(self.multivariateBasis, maxNumPolys, mostSignificant, significanceFactor, True)
    else:
      ErrorMsg('Please provide a valid "strategy" argument. Valid arguments are "Fixed", "Sequential" and "Cleaning"', fileName)

    # Check basis and ask permission to proceed
    if self.checkBasis:
      enum = self.truncatureBasisStrategy.getBasis().getEnumerateFunction()
      for i in range(maxNumPolys):
        print(enum(i))
      print("nPoly = " + str(maxNumPolys))
      inp = raw_input("Do you want to proceed ? (y/n) ")
      if inp in ['Y', 'y', 'Yes', 'yes']:
        pass
      else:
        prt('Exiting', 'red', True)
        sys.exit()
    else:
      pass


  def DefineEvalStrat(self, evalStrategy="Regression", validStrategy="None", nFold=10):
    self.savedItems['CoeffEvaluationStrategy'] = evalStrategy
    self.savedItems['CoeffValidationStrategy'] = validStrategy
    if evalStrategy == "Regression":
      if validStrategy == "None":
        # Default method: The coefficients of the polynomials are found with a standard least-squares algorithm, without validation
        self.coeffEvaluationStrategy = ot.LeastSquaresStrategy()
      else:
        # This algorithm generates a sparse basis from the truncated basis defined above. 
        basisSequenceFactory = ot.LARS()
        if validStrategy == "LOO":
          # Estimates the empirical error on each sub-basis using Leave One Out strategy
          fittingAlgorithm = ot.CorrectedLeaveOneOut()
          # And finally fits the coefficients with a least square algorithm
          approximationAlgorithm = ot.LeastSquaresMetaModelSelectionFactory(basisSequenceFactory, fittingAlgorithm)
          self.coeffEvaluationStrategy = ot.LeastSquaresStrategy(approximationAlgorithm)
        elif validStrategy == "KFold":
          # Estimates the empirical error on each sub-basis using K-fold cross-validation
          # nFold is the number of sets used by the algorithm. Default is 10
          fittingAlgorithm = ot.KFold(nFold)
          # And finally fits the coefficients with a least square algorithm
          approximationAlgorithm = ot.LeastSquaresMetaModelSelectionFactory(basisSequenceFactory, fittingAlgorithm)
          self.coeffEvaluationStrategy = ot.LeastSquaresStrategy(approximationAlgorithm)
        else:
          ErrorMsg('PCE.py/DefineEvalStrat: Please provide a valid "validStrategy" argument. Valid arguments are "LOO", "KFold" and "None"', fileName)
    elif evalStrategy == "Cubature":
      # Cubature-based projection to estimate the coefficients
      self.coeffEvaluationStrategy = ot.IntegrationStrategy()
    else:
      ErrorMsg('Please provide a valid "evalStrategy" argument. Valid arguments are "Regression" and "Cubature"', fileName)


  def SetDistribution(self, distribution):
    self.distribution = distribution


  def _InitWithBlock(self, BLOCK):
    prt('', 'green', self.verbosity)
    prt('Initializing the PCE from the input file.', 'green', self.verbosity)
    
    self.distribution = self.parent.distribution
    self.dimension = self.distribution.getDimension()
    # ----- General information -----
    subBlockNames   = BLOCK.GetSubBlockNames()
    self.ID         = BLOCK.GetCharacterFromKey('ID')
    self.verbosity  = BLOCK.GetCharacterFromKey('outputLevel')
    
    # ----- Polynomial basis definition -----
    if "Polynomial_Basis" in subBlockNames:
      BLOCK.GetSubBlock("Polynomial_Basis")
      keys = BLOCK.GetKeysForBlock()
      if "check" in keys:
        self.checkBasis = BLOCK.GetLogicalFromKey("check")
      else:
        self.checkBasis = False
      if "quasiNorm" in keys:
        quasiNorm = BLOCK.GetRealFromKey("quasiNorm")
      else:
        quasiNorm = 1.0
      if "anisotropyWeights" in keys:
        weights   = BLOCK.GetRealFromKey("anisotropyWeights")
      else:
        weights   = [1.0] * self.dimension #Warning : This imposes that a valid parent experiment has been provided
    else:
      quasiNorm = 1.0
      weights   = [1.0] * self.dimension
    self.DefinePolyBasis(quasiNorm=quasiNorm, weights=weights)

    # ----- Truncature strategy definition -----
    BLOCK.GetSubBlock("Truncature_Strategy")
    keys = BLOCK.GetKeysForBlock()
    strat       = BLOCK.GetCharacterFromKey('strategy')
    if 'maxNumPolys' in keys:
      maxNumPolys = BLOCK.GetIntegerFromKey('maxNumPolys')
      if 'maxTotDegree' in keys:
        WarningMsg('In the presence of the keyword "maxNumPolys", "maxTotDegree" will be ignored', fileName, self.verbosity)
    elif 'maxTotDegree' in keys:
      maxTotDeg = BLOCK.GetIntegerFromKey('maxTotDegree')
      maxNumPolys = self.multivariateBasis.getEnumerateFunction().getStrataCumulatedCardinal(maxTotDeg)
    else:
      ErrorMsg('Either keyword "maxNumPolys" or "maxTotDegree" must be provided', fileName)
    if strat == 'Fixed' or strat == 'Sequential':
      self.DefineTruncatureStrat(strategy=strat, maxNumPolys=maxNumPolys)
    elif strat == 'Cleaning':
      mostSig   = BLOCK.GetIntegerFromKey('mostSignificant')
      signFac   = BLOCK.GetRealFromKey('significanceFactor')
      self.DefineTruncatureStrat(strategy=strat, maxNumPolys=maxNumPolys, mostSignificant=mostSig, significanceFactor=signFac)
    else:
      ErrorMsg('Please provide a valid "strategy" argument. Valid arguments are "Fixed", "Sequential" and "Cleaning"', fileName)

    # ----- Evaluation strategy definition -----
    BLOCK.GetSubBlock("Evaluation_Strategy")
    evalStrat   = BLOCK.GetCharacterFromKey('evalStrategy')
    if evalStrat == "Regression":
      validStrat  = BLOCK.GetCharacterFromKey('validStrategy')
      if validStrat == 'KFold':
        nFold     = BLOCK.GetIntegerFromKey('nFold')
        self.DefineEvalStrat(evalStrategy=evalStrat, validStrategy=validStrat, nFold=nFold)
      else:
        self.DefineEvalStrat(evalStrategy=evalStrat, validStrategy=validStrat)
    else:
      self.DefineEvalStrat(evalStrategy=evalStrat)
