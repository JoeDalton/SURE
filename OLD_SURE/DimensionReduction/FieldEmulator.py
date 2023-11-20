from Utils.I_O import DictFromH5, ErrorMsg, WarningMsg
from PCA import *
from ActiveDirection import *
from SurrogateModelling.PCE import *
from SurrogateModelling.MyKriging import *
from SurrogateModelling.Kriging import *
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import os
from matplotlib.colors import ListedColormap
fileName = 'DimensionReduction/FieldEmulator.py'


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
    #xReSample must be of the shape (nReSample, nUV) # Standard shape of the samplers in SURE
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


  def Build_Single_With_Weights(self, weights):
    reconstruc = self.pca.ComputeTotalFromReduced(weights) #Longest execution time in the routine by a factor 1000

    # Split the meta vector and multiply by the proper scale for each field
    xReconstruc = np.split(reconstruc, self.nField, axis=0)
    for i in range(self.nField):
      xReconstruc[i] *= self.xScale[i]

    # Build back the scaled meta vector and de-center
    reconstruc = np.concatenate(xReconstruc, axis=0)
    reconstruc += self.meanVect[:]
    return reconstruc


  def Evaluate_Single(self, sample):
    # Evaluate surrogates for each mode and build the meta vector
    compressed = np.zeros(self.pca.reducedDim)
    for k in range(self.pca.reducedDim):
      #compressed[k] = self.xSur[k].Evaluate([self.xAD[k].ComputeReducedFromTotal(sample)])[0]
      compressed[k] = self.xSur[k](sample) #DEBUG
    reconstruc = self.pca.ComputeTotalFromReduced(compressed) #Longest execution time in the routine by a factor 1000

    # Split the meta vector and multiply by the proper scale for each field
    xReconstruc = np.split(reconstruc, self.nField, axis=0)
    for i in range(self.nField):
      xReconstruc[i] *= self.xScale[i]

    # Build back the scaled meta vector and de-center
    reconstruc = np.concatenate(xReconstruc, axis=0)
    reconstruc += self.meanVect[:]

    return reconstruc


  def Evaluate_Single_withAD(self, sample):
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


  def PlotField(self, evaluation, fieldID, resr, resx, rlim, xlim, path=None):
    rstep   = 2
    xstep   = 10
    cm = 1.0/2.54  # centimeters in inches
    figsize = (14*cm,10*cm)
    D = 0.00457
    xticks = np.arange(rlim[0], rlim[1], rstep * D)
    xticklabels = np.round(xticks / D).astype(int)
    yticks = np.arange(xlim[0], xlim[1], xstep * D)
    yticklabels = np.round(yticks / D).astype(int)
    fig, ax = plt.subplots(1, 1)#, gridspec_kw={'width_ratios': [1] + [1]*nMode})
    fig.figsize=figsize

    import Utils.PrettyPlots as pp

    def forceAspect(ax,aspect):
      im = ax.get_images()
      extent =  im[0].get_extent()
      ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

    def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
      """Add a vertical color bar to an image plot."""
      divider = axes_grid1.make_axes_locatable(im.axes)
      width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
      pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
      current_ax = plt.gca()
      cax = divider.append_axes("right", size=width, pad=pad)
      plt.sca(current_ax)
      return im.axes.figure.colorbar(im, cax=cax, **kwargs)


    xq = np.linspace(xlim[0],xlim[1],resx)
    rq = np.linspace(rlim[0],rlim[1],resr)
    xcLab = np.zeros(self.nField)

    # Divide the evaluation vector
    xField = np.split(evaluation, self.nField, axis=0) 

    # Plot mean
    j = fieldID
    Q     = np.reshape(xField[j][:],(resr,resx))
    norm  = matplotlib.colors.Normalize(vmin=Q.min(), vmax=Q.max())



    norm = matplotlib.colors.Normalize(vmin=0, vmax=0.001) # for OH
    #Custom colormap to match Cabra 2002 for OH
    cmap = ListedColormap([(0.15,0.18,0.35),(0.15,0.2,0.38),(0.2,0.36,0.58),(0.23,0.42,0.64),(0.25,0.51,0.65),(0.21,0.4,0.24),(0.25,0.49,0.28),(0.54,0.71,0.33),(0.97,0.91,0.37),(0.73,0.33,0.23),(0.67,0.22,0.2),(0.72,0.18,0.2)])
    #TODO: Color as a function of the field name
    #cmap  = 'hot'




    img  = ax.imshow(np.rot90(Q, 1), interpolation='bilinear', cmap=cmap, 
            norm=norm, extent=[rq.min(), rq.max(), xq.min(), xq.max()])
    forceAspect(ax, 0.25)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('r/D [-]')
    ax.set_ylabel('x/D [-]')
    pp.CorrectAxSize(ax, 20)



    cbar = plt.colorbar(img)
    if path is None:
      plt.show()
    else:
      plt.savefig(path, bbox_inches="tight",pad_inches = 0)
    #quit()
    #cbar = plt.colorbar(img, ax=ax[j,1])
    #cbar.set_label(xcLab[j], size=16)
    #cbar.ax[j,1].tick_params(labelsize=12)

    # Print colorbar


    #plt.colorbar(img,fraction=0.046, pad=0.04)


    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax[0])
    #cax = divider.append_axes('right', size='20%', pad=0.05)
    #fig.colorbar(img, cax=cax, orientation='vertical')

    

    #add_colorbar(img)


  def PlotMode(self, fieldID, xVar, resr, resx, rlim, xlim, path=None):

    rstep   = 2
    xstep   = 10
    cm = 1.0/2.54  # centimeters in inches
    #figsize = (14*cm,10*cm)
    figsize = (14*cm,400*cm)
    D = 0.00457
    xticks = np.arange(rlim[0], rlim[1], rstep * D)
    xticklabels = np.round(xticks / D).astype(int)
    yticks = np.arange(xlim[0], xlim[1], xstep * D)
    yticklabels = np.round(yticks / D).astype(int)

    import Utils.PrettyPlots as pp
    nMode = self.pca.reducedDim
    #fig, ax = plt.subplots(1, nMode+2)
    fig, ax = plt.subplots(1, nMode+1) #, gridspec_kw={'width_ratios': [1] + [1]*nMode})
    fig.figsize=figsize

    def forceAspect(ax,aspect):
      im = ax.get_images()
      extent =  im[0].get_extent()
      ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


    xq = np.linspace(xlim[0],xlim[1],resx)
    rq = np.linspace(rlim[0],rlim[1],resr)
    xcLab = np.zeros(self.nField)

    # Divide the mean vector
    xMean = np.split(self.meanVect, self.nField, axis=0) 

    # Plot mean
    j = fieldID
    Q     = np.reshape(xMean[j][:],(resr,resx))
    if xVar[fieldID] == 'OH':
      norm = matplotlib.colors.Normalize(vmin=0, vmax=0.001)
      cmap = ListedColormap([(0.15,0.18,0.35),(0.15,0.2,0.38),(0.2,0.36,0.58),(0.23,0.42,0.64),(0.25,0.51,0.65),(0.21,0.4,0.24),(0.25,0.49,0.28),(0.54,0.71,0.33),(0.97,0.91,0.37),(0.73,0.33,0.23),(0.67,0.22,0.2),(0.72,0.18,0.2)])
    elif xVar[fieldID] == 'T':
      norm  = matplotlib.colors.Normalize(vmin=Q.min(), vmax=Q.max())
      baseCmap=matplotlib.cm.Spectral.reversed()
      midpoint = (1045.0-300.0) / (1500.0-300.0)
      cmap = pp.shiftedColorMap(baseCmap, 0.0, midpoint, 1.0)
    elif xVar[fieldID] == 'z_out':
      cmap  = 'Greys'
      norm  = matplotlib.colors.Normalize(vmin=Q.min(), vmax=Q.max())
    else:
      cmap  = 'Blues'
      norm  = matplotlib.colors.Normalize(vmin=Q.min(), vmax=Q.max())



    img  = ax[0].imshow(np.rot90(Q, 1), interpolation='bilinear', cmap=cmap, 
            norm=norm, extent=[rq.min(), rq.max(), xq.min(), xq.max()])
    forceAspect(ax[0], 0.25)
    ax[0].title.set_text(r'Mean')
    ax[0].set_xticks(xticks)
    ax[0].set_xticklabels(xticklabels)
    ax[0].set_yticks(yticks)
    ax[0].set_yticklabels(yticklabels)
    ax[0].set_xlabel('r/D')
    ax[0].set_ylabel('x/D')
    pp.CorrectAxSize(ax[0], 20)
    #cbar = plt.colorbar(img, ax=ax, orientation='horizontal')
    #cbar.set_label('toto', size=20)

    # Print colorbar


    #plt.colorbar(img,fraction=0.046, pad=0.04)


    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax[0])
    #cax = divider.append_axes('right', size='20%', pad=0.05)
    #fig.colorbar(img, cax=cax, orientation='vertical')

    from mpl_toolkits import axes_grid1
    
    def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
      """Add a vertical color bar to an image plot."""
      divider = axes_grid1.make_axes_locatable(im.axes)
      width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
      pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
      current_ax = plt.gca()
      cax = divider.append_axes("right", size=width, pad=pad)
      plt.sca(current_ax)
      return im.axes.figure.colorbar(im, cax=cax, **kwargs)

    #add_colorbar(img)


    # Plot modes
    xxMode = self.pca.xMode
    for i in range(nMode):
      modeVect = xxMode[i,:] # Here, each block in modeVect is a different field, same as meanVect in the previous paragraph
      
      # Divide the mode vector
      xMode = np.split(modeVect, self.nField, axis=0) 

      # Plot mode
      j = fieldID
      Q     = np.reshape(xMode[j][:],(resr,resx))
      mimax = max(abs(Q.min()), abs(Q.max()))
      Q     = Q/mimax
      #norm  = matplotlib.colors.Normalize(vmin=-mimax, vmax=mimax)
      norm  = matplotlib.colors.Normalize(vmin=-1, vmax=1)
      #norm  = matplotlib.colors.Normalize(vmin=xxMode[:,j].min(), vmax=xxMode[:,j].max())
      
      cmap  = 'RdGy'

      index = i+1
      img   = ax[index].imshow(np.rot90(Q, 1), interpolation='bilinear', cmap=cmap, 
              norm=norm, extent=[rq.min(), rq.max(), xq.min(), xq.max()])
      ax[index].title.set_text(r'Mode ' + str(i+1))
      forceAspect(ax[index], 0.25)
      pp.CorrectAxSize(ax[index], 16)
      ax[0].get_shared_y_axes().join(ax[0], ax[index])
      ax[index].set_yticks([])
      #cbar = plt.colorbar(img, ax=ax[j,1])
      #cbar.set_label(xcLab[j], size=16)
      #cbar.ax[j,1].tick_params(labelsize=12)
  
      ax[index].set_xticks(xticks)
      ax[index].set_xticklabels(xticklabels)
      ax[index].set_xlabel('r/D')
      pp.CorrectAxSize(ax[index], 20)

    fig.tight_layout=True

    if path is None:
      plt.show()
    else:
      plt.savefig(path, bbox_inches='tight')


  def Build_Without_Surrogates(self, xTrainingSample, xTrainingVect):
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

    prt('Field emulator assembled.', 'green', self.verbosity)




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
    self.xSur       = []
    """
    #################################################### With PCE
    for k in range(self.pca.reducedDim):
      self.xSur.append(PCE(checkBasis=False, ID='PCE'+str(k), verbosity=self.verbosity-1, distribution=ot.ComposedDistribution([ot.Normal(0.0,1.0), ot.Uniform(0.0,1.0), ot.Uniform(0.0,1.0)]))) # TODO: Warning, this is for the 3D case only
      self.xSur[k].DefinePolyBasis(quasiNorm=0.5, weights=None)
      self.xSur[k].DefineTruncatureStrat(strategy="Fixed", maxTotDegree=5)
      self.xSur[k].DefineEvalStrat(evalStrategy="Regression", validStrategy="KFold")
      #samples = np.zeros((self.nTrainingSample,1))
      #samples[:, 0] = xxRedSample[k]
      #self.xSur[k].Build(samples, xRedData[k,:])
      self.xSur[k].Build(xTrainingSample, xRedData[k,:])


    """
    #################################################### With myKriging
    for k in range(self.pca.reducedDim):
      # Best cases found with MLE approach
      #if k == 0:
      #  self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], givenTheta=[1.53739247e+00, 1.99411104e-01, 1.00000000e-05], trend='linear'))
      self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], givenTheta=[6.0, 0.5, 0.5], trend='linear'))
      #elif k == 1:
      #  self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], givenTheta=[2.24530129, 0.20332085, 0.00473898], trend='linear'))
      #  #self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], givenTheta=[6.0, 0.4, 0.8], trend='linear'))
      #elif k ==2:
      #  self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], givenTheta=[1.70534681e+00, 1.54061503e-01, 6.70374627e-05], trend='linear'))
      #  #self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], givenTheta=[6.0, 0.4, 0.8], trend='linear'))
      #else:
      #  print("WTF ?? Index too high. I'm reverting to default optimisation")
      #  self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], trend='linear', optistrat='default', optiTarget='MLE'))
        
      

      # Research
      #self.xSur.append(OK_model(ID='K'+str(k), nugget=0.05, af='Matern52', x=xTrainingSample, y=xRedData[k,:], trend='linear', optistrat='default', optiTarget='MLE'))
      #factor = np.array([1.0, 1.0, 5.0])
      #self.xSur[k].theta = np.multiply(self.xSur[k].theta, factor)
      #self.xSur[k].theta[2] = 0.5
      print(self.xSur[k].theta)
      #self.xSur[k].Smoothen(50)


    ##################################################### With OT Kriging
    #for k in range(self.pca.reducedDim):
    #  self.xSur.append(Kriging(ID='K'+str(k), verbosity=self.verbosity-1, distribution=ot.ComposedDistribution([ot.Normal(0.0,1.0), ot.Uniform(0.0,1.0), ot.Uniform(0.0,1.0)]))) # TODO: Warning, this is for the 3D case only
    #  self.xSur[k].DefineTrend('Linear')
    #  self.xSur[k].DefineCovariance('Matern')
    #  self.xSur[k].DefineNoise(1e-9)
    #  self.xSur[k].Build(xTrainingSample, xRedData[k,:])

    prt('Field emulator assembled.', 'green', self.verbosity)


  def Build_WithAD(self, xTrainingSample, xTrainingVect): #Works only for 2D uncertain space at the moment
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
      #self.xAD[k].Build(xTrainingSample, xRedData[k,:]) # Garbage out when xRedData is not monotonous => Mostly garbage for second and higher modes
      self.xAD[k].Build(xTrainingSample, xRedData[0,:]) # Seems to kinda work in the majority of cases, although not rigorous...
      #self._SetArbitraryAD(self.xAD[k], self.ID, k) # TODO: Works only for 2D uncertainty and hand-picked coefficients for known fields...

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
      ErrorMsg(message, fileName)



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
