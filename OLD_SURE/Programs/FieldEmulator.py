from Utils.I_O import DictFromH5, ErrorMsg, WarningMsg
from DimensionReduction.PCA import *
from DimensionReduction.ActiveDirection import *
from SurrogateModelling.PCE import *
import numpy as np
import math

fileName = 'FieldEmulator.py'


class MultipleFieldEmulator():
  # Emulates fields (in the form of vectors) using PCA and surrogate modelling

  ######################
  #     Properties     #
  ######################
  #ID:              ID
  #verbosity:       Verbosity
  #epsilon:         Portion of variance retained by the PCA
  #nTrainingSample: Number of samples in the training set
  #nFeature:        Number of features in the field
  #nField:          Number of fields un the vectors
  #meanVect:        Point-by-point mean vectors of the training fields
  #xScale:          Standard deviations of the centered training fields (one scale by field)
  #pca:             PCA of the training fields
  #xAD:             Collection of active directions to predict the PCA weights
  #xSur:            Collection of surrogates to predict the PCA weights


  ###################
  #     Methods     #
  ###################
  def __init__(self, ID='myFieldEmulator', verbosity=1, epsilon=0.01, nField=1):
    self.ID         = ID
    self.verbosity  = verbosity
    self.epsilon    = epsilon
    self.savedItems = {}
    self.nField     = nField
    self.objectType = 'FieldEmulator'


  def Evaluate_Group(self, xReSample):
    #xTrainingVect must be of the shape (nReSample, nUV) # Standard shape of the samplers in SURE
    nReSample   = xReSample.shape[0]
    xReconstruc = np.zeros((nReSample, self.nFeature*self.nField))
    prt('Reconstructing samples...', 'yellow', self.verbosity)
    pb = ProgBar('Please wait...', maxNumber=nReSample)
    for i in range(nReSample):
      sample = xReSample[i,:]
      xReconstruc[i,:] = self.Evaluate_Single(sample)
      pb.Update()
    pb.Terminate()
    return xReconstruc


  def Evaluate_Single(self, sample):
    # Evaluate surrogates for each mode and build the meta vector
    compressed = np.zeros(self.pca.reducedDim)
    for k in range(self.pca.reducedDim):
      compressed[k] = self.xSur[k].Evaluate([self.xAD[k].ComputeReducedFromTotal(sample)])[0]
    reconstruc = self.pca.ComputeTotalFromReduced(compressed)

    # Split the meta vector and multiply by the proper scale for each field
    xReconstruc = np.split(reconstruc, self.nField, axis=0)
    for i in range(self.nField):
      xReconstruc[i] *= self.xScale[i]

    # Build back the scaled meta vector and de-center
    reconstruc = np.concatenate(xReconstruc, axis=0)
    reconstruc += self.meanVect[:]
  
    # Enjoy !
    return reconstruc      


  def Build(self, xTrainingSample, xTrainingVect):  
    #xTrainingVect must be of the shape (nFeature*nField, nSample)
    prt('Scaling training data...', 'yellow', self.verbosity)
    self.nTrainingSample  = xTrainingVect.shape[1]
    self.nFeature         = xTrainingVect.shape[0]/float(self.nField)
    if self.nFeature.is_integer():
      self.nFeature = int(self.nFeature)
    else:
      ErrorMsg("Inconsistent size of the training set with regard to the number of fields declared. Please try again with consistent input.")
    self.meanVect         = np.mean(xTrainingVect, axis=1)
    centered              = xTrainingVect - self.meanVect[:, np.newaxis]
    xCentered             = np.split(centered, self.nField, axis=0)
    self.xScale           = np.zeros(self.nField)
    xScaled               = []
    for i in range(self.nField):
      self.xScale[i]      = np.std(xCentered[i])
      xScaled.append(xCentered[i] / self.xScale[i])
    scaled                = np.concatenate(xScaled, axis=0)

    prt('Computing PCA modes...', 'yellow', self.verbosity)
    self.pca        = PCA('myPCA', scaled, epsilon=self.epsilon, parent=self, verbosity=self.verbosity-1)

    prt('Constructing reduced dataset...', 'yellow', self.verbosity)
    xRedData        = np.zeros((self.pca.reducedDim, self.nTrainingSample))
    for i in range(self.nTrainingSample):
      xRedData[:,i] = self.pca.ComputeReducedFromTotal(scaled[:,i])

    prt('Constructing surrogates...', 'yellow', self.verbosity)
    self.xAD        = []
    xxRedSample     = []
    self.xSur       = []

    for k in range(self.pca.reducedDim):
      #TODO: Make this a lot more modular, with a choice for Active direction or not, type of surrogate...
      self.xAD.append(ActiveDirection(ID='AD'+str(k), totalDim=xTrainingSample.shape[1], verbosity=1))
      #xAD[k].Build(xTrainingSample, xRedData[k,:]) # Garbage out when xRedData is not monotonous => Mostly garbage for second and higher modes
      #xAD[k].Build(xTrainingSample, xRedData[0,:]) # Seems to kinda work in the majority of cases, although not rigorous...
      self._SetArbitraryAD(self.xAD[k], self.ID, k) # TODO: Works only for 2D uncertainty and hand-picked coefficients for known fields...

      xxRedSample.append(self.xAD[k].ComputeReducedFromTotal(xTrainingSample))
      #self.xSur.append(PCE(ID='PCE'+str(k), verbosity=self.verbosity-1, distribution=ot.ComposedDistribution([ot.Normal(0.0,1.0)]))) # TODO: Not Modular at all
      self.xSur.append(PCE(ID='PCE'+str(k), verbosity=self.verbosity-1, distribution=ot.ComposedDistribution([ot.Uniform(-1.0,0.0)]))) # TODO: Not Modular at all
      self.xSur[k].DefinePolyBasis(quasiNorm=1.0, weights=None) 
      self.xSur[k].DefineTruncatureStrat(strategy="Fixed", maxTotDegree=5)
      self.xSur[k].DefineEvalStrat(evalStrategy="Regression", validStrategy="KFold")
      samples = np.zeros((self.nTrainingSample,1))
      samples[:, 0] = xxRedSample[k]
      self.xSur[k].Build(samples, xRedData[k,:])

    prt('Field emulator assembled.', 'green', self.verbosity)


  def _SetArbitraryAD(self, AD, variable, k):
    if variable == 'OH':
      if k == 0:
        angleDeg  = 1
      elif k==1:
        angleDeg  = 3
      elif k==2:
        angleDeg  = 2
      else:
        message = 'Unknown case, we set a default case like the first AD'
        WarningMsg(message, fileName, self.verbosity)
        angleDeg = 1
      angle     = angleDeg*math.pi/180.0
      AD.xWeight = [-math.sin(angle), -math.cos(angle)]
      
    else: 
      message = 'Unknown variable. Exiting...'
      ErrorMsg(message, fileName, self.verbosity)



class FieldEmulator():
  # Emulates fields (in the form of vectors) using PCA and surrogate modelling

  ######################
  #     Properties     #
  ######################
  #ID:              ID
  #verbosity:       Verbosity
  #epsilon:         Portion of variance retained by the PCA
  #nTrainingSample: Number of samples in the training set
  #nFeature:        Number of features in the field
  #meanVect:        Point-by-point mean of the training fields
  #scale:           Standard deviation of the centered training fields
  #pca:             PCA of the training fields
  #xAD:             Collection of active directions to predict the PCA weights
  #xSur:            Collection of surrogates to predict the PCA weights


  ###################
  #     Methods     #
  ###################
  def __init__(self, ID='myFieldEmulator', verbosity=1, epsilon=0.01):
    self.ID         = ID
    self.verbosity  = verbosity
    self.epsilon    = epsilon
    self.savedItems = {}
    self.objectType = 'FieldEmulator'


  def Evaluate_Group(self, xReSample):
    #xTrainingVect must be of the shape (nReSample, nUV) # Standard shape of the samplers in SURE
    nReSample   = xReSample.shape[0]
    xReconstruc = np.zeros((nReSample, self.nFeature))
    prt('Reconstructing samples...', 'yellow', self.verbosity)
    pb = ProgBar('Please wait...', maxNumber=nReSample)
    for i in range(nReSample):
      sample = xReSample[i,:]
      xReconstruc[i,:] = self.Evaluate_Single(sample)
      pb.Update()
    pb.Terminate()
    return xReconstruc


  def Evaluate_Single(self, sample):
    compressed = np.zeros(self.pca.reducedDim)
    for k in range(self.pca.reducedDim):
      compressed[k] = self.xSur[k].Evaluate([self.xAD[k].ComputeReducedFromTotal(sample)])[0]
    reconstruc = self.pca.ComputeTotalFromReduced(compressed)
    reconstruc *= self.scale
    reconstruc += self.meanVect[:]
    return reconstruc      


  def Build(self, xTrainingSample, xTrainingVect):  
    #xTrainingVect must be of the shape (nFeature, nSample)
    prt('Scaling training data...', 'yellow', self.verbosity)
    self.nTrainingSample  = xTrainingVect.shape[1]
    self.nFeature         = xTrainingVect.shape[0]
    self.meanVect         = np.mean(xTrainingVect, axis=1)
    centered              = xTrainingVect - self.meanVect[:, np.newaxis]
    self.scale            = np.std(centered)
    scaled                = centered / self.scale

    prt('Computing PCA modes...', 'yellow', self.verbosity)
    self.pca        = PCA('myPCA', scaled, epsilon=self.epsilon, parent=self, verbosity=self.verbosity-1)

    prt('Constructing reduced dataset...', 'yellow', self.verbosity)
    xRedData        = np.zeros((self.pca.reducedDim, self.nTrainingSample))
    for i in range(self.nTrainingSample):
      xRedData[:,i] = self.pca.ComputeReducedFromTotal(scaled[:,i])

    prt('Constructing surrogates...', 'yellow', self.verbosity)
    self.xAD        = []
    xxRedSample     = []
    self.xSur       = []

    for k in range(self.pca.reducedDim):
      #TODO: Make this a lot more modular, with a choice for Active direction or not, type of surrogate...
      self.xAD.append(ActiveDirection(ID='AD'+str(k), totalDim=xTrainingSample.shape[1], verbosity=1))
      #xAD[k].Build(xTrainingSample, xRedData[k,:]) # Garbage out when xRedData is not monotonous => Mostly garbage for second and higher modes
      #xAD[k].Build(xTrainingSample, xRedData[0,:]) # Seems to kinda work in the majority of cases, although not rigorous...
      self._SetArbitraryAD(self.xAD[k], self.ID, k) # TODO: Works only for 2D uncertainty and hand-picked coefficients for known fields...

      xxRedSample.append(self.xAD[k].ComputeReducedFromTotal(xTrainingSample))
      #self.xSur.append(PCE(ID='PCE'+str(k), verbosity=self.verbosity-1, distribution=ot.ComposedDistribution([ot.Normal(0.0,1.0)]))) # TODO: Not Modular at all
      self.xSur.append(PCE(ID='PCE'+str(k), verbosity=self.verbosity-1, distribution=ot.ComposedDistribution([ot.Uniform(-1.0,0.0)]))) # TODO: Not Modular at all
      self.xSur[k].DefinePolyBasis(quasiNorm=1.0, weights=None) 
      self.xSur[k].DefineTruncatureStrat(strategy="Fixed", maxTotDegree=5)
      self.xSur[k].DefineEvalStrat(evalStrategy="Regression", validStrategy="KFold")
      samples = np.zeros((self.nTrainingSample,1))
      samples[:, 0] = xxRedSample[k]
      self.xSur[k].Build(samples, xRedData[k,:])

    prt('Field emulator assembled.', 'green', self.verbosity)


  def _SetArbitraryAD(self, AD, variable, k):
    if variable == 'OH':
      if k == 0:
        angleDeg  = 1
      elif k==1:
        angleDeg  = 3
      elif k==2:
        angleDeg  = 2
      else:
        message = 'Unknown case, we set a default case like the first AD'
        WarningMsg(message, fileName, self.verbosity)
        angleDeg = 1
      angle     = angleDeg*math.pi/180.0
      AD.xWeight = [-math.sin(angle), -math.cos(angle)]
      
    else: 
      message = 'Unknown variable. Exiting...'
      ErrorMsg(message, fileName, self.verbosity)
