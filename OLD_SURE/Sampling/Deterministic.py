from Utils.I_O import prt, ErrorMsg, WarningMsg
from Numerics.Algebra import norm1
from Numerics.Algebra import CheckVectorsEquality
from Sampling  import Sampling
import openturns as ot
import sys
import math
from scipy.special import erfinv
import numpy as np
#import pickle


class DetSampling(Sampling):

  ######################
  #     Properties     #
  ######################
  #xPoint             = None # np.array Already defined in parent class
  xWeight            = None # np.array
  fileName           = ''
  ###################
  #     Methods     #
  ###################
  def __init__(self, parent=None,verbosity=-1):
    self.parent     = parent
    self.fileName = 'Sampling/Deterministic.py/DetSampling'
    self.verbosity  = verbosity
 

  #def ToOT(self): #TODO: Remove this ?
  #  sample = ot.Sample(self.xPoint)
  #  result = ot.FixedExperiment(sample) 
  #  # WARNING ! the FixedExperiment cannot be set up with the weights.
  #  # This is because the constructor does not accept negative weights as it interprets them as negative probabilities
  #  # ** sigh **
  #  # The weights are thus automatically set to 1/size(sample) = 1/len(xPoint) = 1/len(xWeight)
  #  # To prepare the proper IntegrationStrategy for the PCE, the result given by the solver for sample[i] 
  #  # must be multiplied by len(xWeight) * xWeight[i]
  #  return result

 
  def __add__(self, other):
    return self._AddSub(other, 1.0)


  def __sub__(self, other):
    return self._AddSub(other, -1.0)


  def _AddSub(self, other, coeff):
    if isinstance(self, CubatureSampling):
      result        = CubatureSampling(self.parent, xRule=self.xRule, xDistributionType=self.xDistributionType, \
                                       xxDistributionParam=self.xxDistributionParam, maxLevel=self.maxLevel, \
                                       xLevel=self.xLevel, isSparse = self.isSparse, verbosity=self.verbosity-1)
    elif isinstance(self, QuadratureSampling):
      result        = QuadratureSampling(self.parent, self.rule, self.transformFunc, \
                                         xDistributionParam=self.xTransformParam, level=self.level, \
                                         nPoint=self.nPoint, verbosity=self.verbosity-1)
    else:
      result        = DetSampling()
    result.xPoint   = self.xPoint
    result.xWeight  = self.xWeight

    if result.xPoint is None:
      if other.xPoint is not None:
        result.xPoint  = np.array(other.xPoint)
        result.xWeight = coeff * np.array(other.xWeight)
      return result
    if other.xPoint is None:
      return result
    
    # For each point of self, add or substract the weights of the corresponding point of other (if they exist)
    for i in range(len(result.xPoint)):
      indices = []
      for k in range(len(other.xPoint)):
        if CheckVectorsEquality(result.xPoint[i], other.xPoint[k], rtol=1e-15):
          indices.append(k)

      if len(indices) == 0:
        pass #Nothing to do in this case.
      elif len(indices) > 1:
        ErrorMsg('"indices" cannot have a length superior than 1. Please check your quadrature methods', fileName)
      else:
        # A point equal to result.xPoint[i] has been found in other.xPoint. We substract the weight of the later to the weight of the former
        k = indices[0]
        result.xWeight[i] += coeff * other.xWeight[k] 
    
    # For each point of other, if they are not in self, add them in result with opposite weight
    for j in range(len(other.xPoint)):
      indices = []
      for k in range(len(result.xPoint)):
        if CheckVectorsEquality(result.xPoint[k], other.xPoint[j], rtol=1e-15):
          indices.append(k)
      if len(indices) == 0:
        result.xPoint  = np.append(result.xPoint, [other.xPoint[j]], axis=0)
        result.xWeight = np.append(result.xWeight, [coeff * other.xWeight[j]], axis=0)
    return result











class CubatureSampling(DetSampling):

  ######################
  #     Properties     #
  ######################
  nPoint              = 0     # Integer
  maxLevel            = 0     # Integer
  xLevel              = []    # List of integers
  xRule               = []    # List of strings
  xDistributionType   = []    # List of strings
  xxDistributionParam = []    # List of lists of integers
  xQuadDrawer         = []    # List of quadrature samplings
  isSparse            = False # Bool

  ###################
  #     Methods     #
  ###################
  def __init__(self, parent, **kwargs):
    # Set parent experiment
    self.parent = parent
    self.fileName = 'Sampling/Deterministic.py/CubatureSampling'

    # Misc
    self.xQuadDrawer = []
    if 'verbosity' in kwargs:
      self.verbosity = kwargs['verbosity']
    else:
      self.verbosity = -1

    # Set variables distribution
    if 'xDistributionType' in kwargs and 'xxDistributionParam' in kwargs:
      self.xDistributionType    = kwargs['xDistributionType']
      self.xxDistributionParam  = kwargs['xxDistributionParam']
      dim = len(self.xDistributionType)
      if len(self.xxDistributionParam) != dim:
        ErrorMsg('"xxDistribution" must have the same length as xDistributionType', self.fileName)
      if 'distribution' in kwargs:
        WarningMsg('Argument "distribution" is ignored in the presence of "xDistributionType" and "xxDistributionParam"', self.fileName, self.verbosity)
    elif 'distribution' in kwargs:
      distribution        = kwargs['distribution']
      dim                 = distribution.getDimension()
      xDistributionType   = []
      xxDistributionParam = []
      for i in range(dim):
        xDistributionType.append(distribution.getMarginal(i).getName())
        xxDistributionParam.append([0.0,1.0])
      self.xDistributionType    = xDistributionType
      self.xxDistributionParam  = xxDistributionParam
    else:
      ErrorMsg('You must provide either one of the two following groups of arguments: "distribution" or ("xDistributionType" and "xxDistributionParam")', self.fileName)

    # Set cubature rules
    if 'xRule' in kwargs:
      self.xRule = kwargs['xRule']
      if len(self.xRule) != dim:
        ErrorMsg('"xRule" must have the same dimension as the distribution', self.fileName)
      if 'rule' in kwargs:
        WarningMsg('Argument "rule" is ignored in the presence of "xRule"', self.fileName, self.verbosity)
    elif 'rule' in kwargs:
      self.xRule = [kwargs['rule']] * dim
    else:
      ErrorMsg('You must provide either one of the two following arguments: "rule" or "xRule"', self.fileName)

    # Set cubature levels
    if 'xLevel' in kwargs:
      self.xLevel = kwargs['xLevel']
      if isinstance(self.xLevel, int):
        self.xLevel = [self.xLevel]
      self.maxLevel = max(self.xLevel)
      if len(self.xRule) != dim:
        ErrorMsg('"xLevel" must have the same dimension as the distribution', self.fileName)
      if 'level' in kwargs:
        WarningMsg('Argument "level" is ignored in the presence of "xLevel"', self.fileName, self.verbosity)
    elif 'level' in kwargs:
      self.xLevel   = [kwargs['level']] * dim
      self.maxLevel = kwargs['level']
    elif 'nPointPerDim' in kwargs:
      self.nPoint   = kwargs['nPointPerDim'] ** dim
    else:
      ErrorMsg('You must provide either one of the three following arguments: "nPointPerDim" (for Regular sampling) "level" or "xLevel"', self.fileName)

    # Set sparsity
    if 'isSparse' in kwargs:
      self.isSparse = kwargs['isSparse']
    else:
      WarningMsg('Default: No sparsity', self.fileName, self.verbosity)
      self.isSparse = False


  def Draw(self):
    if self.xRule[0] == 'Regular':
      self._DrawRegular()
    else:
      if self.isSparse:
        self._DrawSmolyak()
      else:
        self._DrawFull()


  def _DrawSmolyak(self):
    xMultiIndex = self._ComputeAdmissibleMultiIndices()
    for multiIndex in xMultiIndex:
      temp = self + self._ComputeTensDeltaL(multiIndex)
      self.xPoint = temp.xPoint
      self.xWeight = temp.xWeight


  def _ComputeAdmissibleMultiIndices(self):
    dim = len(self.xLevel)

    # Creating the precursor lists of levels for each coordinates
    xxCoorLevel=[]
    for l in self.xLevel:
      xxCoorLevel.append(range(0, l+1))

    # Enumerating all possible multiIndices given xLevel
    xMultiIndex = []
    for i in range(len(xxCoorLevel[0])):
      xMultiIndex.append([xxCoorLevel[0][i]])

    # Constructing the vector of multiIndices by adding iteratively the new levels on the whole previous vector
    for i in range(1, dim):
      newXMultiIndex = []
      # Retrieving the list of indices for variable i
      currentXLevel  = xxCoorLevel[i]

      # Core of the procedure:
      for j in range(len(currentXLevel)):
        for k in range(len(xMultiIndex)):
          # Creating new multiIndices including the new dimension
          tempSample = xMultiIndex[k] + [currentXLevel[j]]
          newXMultiIndex.append(tempSample)
      xMultiIndex = list(newXMultiIndex) #copy as value and not as reference

    # Eliminate all multiIndices non-admissible given maxLevel
    tempXMultiIndex = list(xMultiIndex)
    xMultiIndex     = []
    for multiIndex in tempXMultiIndex:
      if norm1(multiIndex)<=self.maxLevel:
        xMultiIndex.append(multiIndex)

    # Return result
    return xMultiIndex


  def _ComputeTensDeltaL(self, multiIndex):
    TensDeltaL = CubatureSampling(self.parent, xRule=self.xRule, xDistributionType=self.xDistributionType, \
                                  xxDistributionParam=self.xxDistributionParam, maxLevel=self.maxLevel, \
                                  xLevel=self.xLevel, isSparse = self.isSparse, verbosity=self.verbosity-1)
    if len(multiIndex) != len(self.xRule):
      ErrorMsg('Argument "multiIndex" must be of same size than "self.xRule"', self.fileName)
    else:
      # Compute the partial quadrature samplings i of level multiIndex[i] 
      xDeltaL = []
      for i in range(len(multiIndex)):
        xDeltaL.append(self._ComputeDeltaL(i, multiIndex[i]))
      # The result is the tensorisation of these quadratures
      TensDeltaL.xPoint, TensDeltaL.xWeight = self._TensorizeQuads(xDeltaL)
    return TensDeltaL


  def _ComputeDeltaL(self, coord, level):
    if level == 0:
      # Delta0 = Q0
      Q0Drawer = QuadratureSampling(self, self.xRule[coord], self.xDistributionType[coord], \
                                    xDistributionParam=self.xxDistributionParam[coord], level=level)
      Q0Drawer.Draw()
      deltaL = Q0Drawer
    else:
      #DeltaL = Ql - Ql-1
      QlDrawer          = QuadratureSampling(self, self.xRule[coord], self.xDistributionType[coord], \
                                             xDistributionParam=self.xxDistributionParam[coord], level=level)
      QlMinusOneDrawer  = QuadratureSampling(self, self.xRule[coord], self.xDistributionType[coord], \
                                             xDistributionParam=self.xxDistributionParam[coord], level=level-1)
      QlDrawer.Draw()
      QlMinusOneDrawer.Draw()
      deltaL = QlDrawer - QlMinusOneDrawer
    return deltaL


  def _DrawFull(self):
    # Initialize quadrature samplings
    self.xQuadDrawer = []
    for i in range(len(self.xRule)):
      quadDrawer = QuadratureSampling(self, self.xRule[i], self.xDistributionType[i], \
                                      xDistributionParam=self.xxDistributionParam[i], level=self.xLevel[i])
      self.xQuadDrawer.append(quadDrawer)
    # Draw quadrature points
    for quad in self.xQuadDrawer:
      quad.Draw()
    # Tensorize the quadratures
    self.xPoint, self.xWeight = self._TensorizeQuads(self.xQuadDrawer)


  def _DrawRegular(self):
    # Initialize quadrature samplings
    self.xQuadDrawer = []
    for i in range(len(self.xRule)):
      quadDrawer = QuadratureSampling(self, self.xRule[i], self.xDistributionType[i], \
                                      xDistributionParam=self.xxDistributionParam[i], nPoint=int(self.nPoint ** (1.0/(len(self.xRule)))))
      self.xQuadDrawer.append(quadDrawer)
    # Draw quadrature points
    for quad in self.xQuadDrawer:
      quad.Draw()
    # Tensorize the quadratures
    self.xPoint, self.xWeight = self._TensorizeQuads(self.xQuadDrawer)


  def _TensorizeQuads(self, xQuad):
    # Initialising of the sample array (list of list for now)
    xSample = []
    xWeight = []
    for i in range(len(xQuad[0].xPoint)):
      xSample.append([xQuad[0].xPoint[i]])    
    xWeight = xQuad[0].xWeight.tolist()
    
    # Constructing the vector of samples by adding iteratively the new coordinate on the whole previous vector
    for i in range(1, len(xQuad)):
      newXSample = []
      newXWeight = []
      # Retrieving the list of points for variable i
      currentXCoord  = xQuad[i].xPoint.tolist()
      currentXWeight = xQuad[i].xWeight.tolist()
      
      # Core of the procedure:
      for j in range(len(currentXCoord)):
        for k in range(len(xSample)):
          # Creating new points with the coordinates of the current variable
          tempSample = xSample[k] + [currentXCoord[j]]
          newXSample.append(tempSample)
          # Creating the weight of the points
          tempWeight = xWeight[k] * currentXWeight[j]
          newXWeight.append(tempWeight)
      xSample = list(newXSample) #copy as value and not as reference
      xWeight = list(newXWeight)
    return np.array(xSample), np.array(xWeight)










class QuadratureSampling(DetSampling):

  ######################
  #     Properties     #
  ######################
  level             = 0     # Integer
  xTransformParam   = []    # List of integers
  transformFunc     = ''    # String
  rule              = ''    # String
  nPoint            = 0     # Integer
  xRUPoint          = None  # np.array

  ###################
  #     Methods     #
  ###################
  def __init__(self, parent, rule, distributionType, xDistributionParam=[0.0,1.0], level=0, nPoint=0, verbosity=-1):
    self.parent       = parent
    self.rule         = rule
    self.level        = level
    self.verbosity    = verbosity
    self.fileName = 'Sampling/Deterministic.py/QuadratureSampling'

    if self.rule in ['SecondFejer', 'Clenshaw-Curtis', 'Simpson']:
      self._ComputeNPoint()
    else:
      if level != 0:
        WarningMsg('For the selected rule, the "level" argument is not relevant. It will be ignored.', self.fileName, self.verbosity)
      if nPoint == 0:
        ErrorMsg('For the selected rule, the "nPoint" argument must be provided.', self.fileName)
      else:
        self.nPoint = nPoint
  
    # Prepare transformation from the U[-1,1] sampling to the real distribution
    if distributionType in ['Uniform', 'LogUniform', 'Normal', 'LogNormal']:
      if len(xDistributionParam) != 2:
        ErrorMsg('For the selected "distributionType", "xDistributionParam" must be of length 2', self.fileName)
    elif distributionType == 'Dirac':
      if len(xDistributionParam) != 1:
        ErrorMsg('For the selected "distributionType", "xDistributionParam" must be of length 1', self.fileName)
    else:
      ErrorMsg('The provided argument for "distributionType" is not supported at the moment. Valid arguments are "Uniform", "Dirac", "Normal" and "LogNormal"', self.fileName)
    self.xTransformParam = xDistributionParam
    self.transformFunc   = distributionType

      


  def _ComputeNPoint(self):
    if self.rule == 'SecondFejer':
      self.nPoint = 2**(self.level+1) - 1
    elif self.rule == 'Clenshaw-Curtis':
      self.nPoint = 2**self.level + 1
    elif self.rule == 'Simpson':
      self.nPoint = 2 * self.level + 1


  def Draw(self):
    # Clear
    self.xWeight  = np.zeros(self.nPoint)
    self.xRUPoint = np.zeros(self.nPoint)
    self.xPoint   = np.zeros(self.nPoint)
    # Draw on the reduced-Uniform interval [-1,1] according to the selected rule and number of points
    if self.rule == 'SecondFejer':
      self.DrawSecondFejer()
    elif self.rule == 'Clenshaw-Curtis':
      self.DrawClenshawCurtis()
    elif self.rule == 'Simpson':
      self.DrawSimpson()
    elif self.rule == 'Regular':
      self.DrawRegular()
    else:
      ErrorMsg('Invalid "rule" argument. Valid arguments are "SecondFejer", "Clenshaw-Curtis", "Simpson" and "Regular"', self.fileName)
    # Transform the points to fit them to their initial distribution
    self._Transform() #TODO: Maybe remove this and put the transformation directly in the program ?. For cubatures, transform only after the tensorization


  def DrawSecondFejer(self, xMin=-1.0, xMax=1.0):
    if self.nPoint == 1:
      self.xPoint[0]  = (xMin + xMax) / 2
      self.xWeight[0] = 1.0
    else:
      N = self.nPoint + 1
      xMin *= 1.0
      xMax *= 1.0
      for i in range(1, self.nPoint+1):
        xi = math.cos(i * math.pi / N)              # Find evaluation point in [-1; 1]
        xi = ((xMax - xMin) * xi + xMax + xMin) / 2 # Convert to [xMin, xMax]
        Ki = 0.0
        for m in range(1, int(math.floor(N / 2)) + 1):
            Ki += math.sin((2 * m - 1) * i * math.pi / N) / (2 * m - 1)
        wi = 4 * math.sin(i * math.pi / N) * Ki / N
        wi *= 1/(xMax - xMin)
        
        self.xWeight[i-1] = wi
        self.xRUPoint[i-1] = xi


  def DrawClenshawCurtis(self, xMin=-1.0, xMax=1.0):
    if self.nPoint == 1:
      self.xPoint[0]  = (xMin + xMax) / 2
      self.xWeight[0] = 1.0
    else:
      N = self.nPoint -1
      xMin *= 1.0
      xMax *= 1.0
      for i in range(0, self.nPoint):
        xi = math.cos(i * math.pi / N)              # Find evaluation point in [-1; 1]
        xi = ((xMax - xMin) * xi + xMax + xMin) / 2 # Convert to [xMin, xMax]
        K = 0.0
        if i == 0 or i == N:
            ci = 1
        else:
            ci = 2
        for j in range(1, int(math.floor(N/2))+1):
            if j == N/2:
                bj = 1
            else:
                bj = 2
            K += bj * math.cos(2*j*i*math.pi/N) / (4*j*j - 1)
        wi = ci * (1-K)/N
        wi *= 1/(xMax-xMin)

        self.xWeight[i] = wi
        self.xRUPoint[i] = xi


  def DrawSimpson(self, xMin=-1.0, xMax=1.0):
    if self.nPoint == 1:
      self.xPoint[0]  = (xMin + xMax) / 2
      self.xWeight[0] = 1.0
    else:
      xMin *= 1.0
      xMax *= 1.0
      step = (xMax-xMin) / (self.nPoint-1)
      for i in range(self.nPoint):
        # This should work
        xi = xMin + i * step
        if i == 0 or i == self.nPoint -1:
          wi = (xMax - xMin) / (self.level * 1.0)
        elif i%2 == 1: 
          wi = 4.0 * (xMax - xMin) / (self.level * 1.0)
        else:
          wi = 2.0 * (xMax - xMin) / (self.level * 1.0)
        self.xWeight[i] = wi / 12.0 # TODO: Maybe check this, but it seems to work just fine
        self.xRUPoint[i] = xi


  def DrawRegular(self, xMin=-1.0, xMax=1.0):
    xMin *= 1.0
    xMax *= 1.0
    step = (xMax-xMin) / (self.nPoint-1)
    for i in range(self.nPoint):
      xi = xMin + i * step
      wi = 1.0/self.nPoint
      self.xWeight[i] = wi
      self.xRUPoint[i] = xi








 
################################# TODO: Maybe move this to DetSampling and make it compatible with scalars and arrays (see test function below). 
################################# Also take transformFunc, xTranformParam, xPoint and xRUPoint as optional parameters. If not provided, take self.*
 
  def _Transform(self):
    if self.transformFunc == 'Dirac':
      self._TransformDirac()
    elif self.transformFunc == 'Uniform':
      self._TransformUniform()
      #elif self.transformFunc == 'LogUniform':
      #  self._TransformLogUniform()
    elif self.transformFunc == 'Normal':
      self._TransformNormal()
    elif self.transformFunc == 'LogNormal':
      self._TransformLogNormal()
    else:
      ErrorMsg('Invalid "self.tranformFunc" argument. Valid arguments are "Dirac", "Uniform", "Normal",  and "LogNormal"', self.fileName)

  def _TransformDirac(self):
    # Simply collapse every result on the dirac value
    value = 1.0 * self.xTransformParam[0]
    self.xPoint[:] = value

  def _TransformUniform(self):
    # Simply translate and stretch the result
    xMin = 1.0 * self.xTransformParam[0]
    xMax = 1.0 * self.xTransformParam[1]
    self.xPoint[:] = xMin + (self.xRUPoint[:]+1) * (xMax-xMin)/2
    self.xWeight[:] = self.xWeight[:] * (xMax-xMin) # TODO: Sure about that ???

  #def _TransformLogUniform(self):
  #  # Simply translate and stretch the result
  #  xMin = 1.0 * self.xTransformParam[0]
  #  xMax = 1.0 * self.xTransformParam[1]
  #  self.xPoint[:] = np.exp(math.log(xMin) + (self.xRUPoint[:]+1) * (math.log(xMax)-math.log(xMin))/2)
  #  #TODO: What about the weights ?

  def _TransformNormal(self):
    # Use the inversion method, with the quantile function of the normal distribution adapted for initial sampling in [-1,1]
    mu    = 1.0 * self.xTransformParam[0]
    sigma = 1.0 * self.xTransformParam[1]
    self.xPoint[:] = mu + sigma * math.sqrt(2) * erfinv(self.xRUPoint[:])

  def _TransformLogNormal(self):
    # Use the inversion method, with the quantile function of the log-normal distribution adapted for initial sampling in [-1,1]
    mu    = 1.0 * self.xTransformParam[0]
    sigma = 1.0 * self.xTransformParam[1]
    self.xPoint[:] = np.exp(mu + sigma * math.sqrt(2) * erfinv(self.xRUPoint[:]))

  #def func_for_scalars_or_vectors(x):
  #  x = np.asarray(x)
  #  scalar_input = False
  #  if x.ndim == 0:
  #    x = x[None]  # Makes x 1D
  #    scalar_input = True
  #  
  #  # The magic happens here
  #  
  #  if scalar_input:
  #    return np.squeeze(ret)
  #  return ret

##################################################
##################################################

