#!/usr/bin/env python2.7

import os
import sys
import h5py
import Link
import argparse
import numpy                                                      as np
import openturns                                                  as ot
import matplotlib.pyplot                                          as plt
#import Utils.PrettyPlots                                          as pp
from Utils.I_O                        import prt
from Utils.I_O                        import DumpObject
from Utils.I_O                        import InputHandler         as IH
from Utils.I_O                        import GetObjectType
from scipy.stats                      import wasserstein_distance as wdist
from Sampling.MC                      import MC
from Numerics.Algebra                 import norm2, norminf
from Experiment.Experiment            import Experiment
from SurrogateModelling.PCE           import *
from SurrogateModelling.OneDSpline    import *
from SurrogateModelling.Kriging       import *
from DimensionReduction.Preprocessing import Preprocessor


import matplotlib as mpl
font = {'family' : 'sans',
        'serif': 'helvet',
        'weight' : 'bold',
        'size'   : 20}

mpl.rc('font', **font)
mpl.rc('text', usetex=True)



#TODO: Better input file with independent blocks for "Data", "Experiment", "Surrogate",
# and such. Links to these blocks in the Validation and reference blocks

# ----- Parsing arguments -----
parser      = argparse.ArgumentParser(description="SurrogateChecker")
parser.add_argument('-i',	          required=False,	type=str,   help = "Relative path of input file")
args        = parser.parse_args()

# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'SurrogateChecker.input'
  try:
    INPUT = IH(inputFile, modify=False)
  except:
    prt('Input File does not exist', 'None', True)
    prt('  --> Copy an example in ' + inputFile, 'None', True)
    os.system("cp " + Link.sureHome + '/InputFiles/' + inputFile + ' .')
    sys.exit()


#TODO: Maybe add the possibility check several surrogates at the same time (useful for projected PCEs of different level/degree for example)

##############################
#     INPUT FILE READING     #
##############################
# ----- Surrogate-related inputs -----
INPUT.GetBlock('Surrogate')
keys            = INPUT.GetKeysForBlock()
surFile         = INPUT.GetCharacterFromKey('SurrogateFile')
surPath         = INPUT.GetCharacterFromKey('SurrogatePath')
surDataIPath    = INPUT.GetCharacterFromKey('InputsPath')
surDataOPath    = INPUT.GetCharacterFromKey('OutputsPath')
ExpeFile        = INPUT.GetCharacterFromKey('ExperimentFile')
ExpePath        = INPUT.GetCharacterFromKey('ExperimentPath')
surType         = GetObjectType(surFile, surPath)
sumFile         = INPUT.GetCharacterFromKey('SummaryPlotFile')
modFile         = ""
if "ModelPlotFile" in keys:
  modFile       = INPUT.GetCharacterFromKey('ModelPlotFile')
if "PlotKDE" in keys:
  surKDE        = INPUT.GetLogicalFromKey('PlotKDE')
else:
  surKDE        = False
if "SmoothenKDE" in keys:
  surKDESmooth  = INPUT.GetLogicalFromKey('SmoothenKDE')
else:
  surKDESmooth  = False
if "PreprocessingSurrogate" in keys:
  isSurPreProc  = True
  surPreProc    = INPUT.GetCharacterFromKey('PreprocessingSurrogate')
else:
  isSurPreProc  = False
if "PreprocessingOriPoints" in keys:
  isOriPreProc  = True
  oriPreProc    = INPUT.GetCharacterFromKey('PreprocessingOriPoints')
else:
  isOriPreProc  = False
if 'SurrogateOnOriPoints' in keys:
  surOnOri = INPUT.GetLogicalFromKey('SurrogateOnOriPoints')
else:
  surOnOri = True

# ----- Check if reference and validation are provided -----
isValidation    = False
isReference     = False
try:
  INPUT.GetBlock('Validation')
  isValidation  = True
  valFile       = INPUT.GetCharacterFromKey('ValidationFile')
  valDataIPath  = INPUT.GetCharacterFromKey('InputsPath')
  valDataOPath  = INPUT.GetCharacterFromKey('OutputsPath')
  keys = INPUT.GetKeysForBlock()
  if "Preprocessing" in keys:
    isValPreProc = True
    valPreProc = INPUT.GetCharacterFromKey('Preprocessing')
  else:
    isValPreProc = False
except:
  pass
try:
  INPUT.GetBlock('Reference')
  isReference   = True
  keys = INPUT.GetKeysForBlock()
  pdfFile       = INPUT.GetCharacterFromKey('PDFPlotFile')
  refFile       = INPUT.GetCharacterFromKey('ReferenceFile')
  refDataIPath  = INPUT.GetCharacterFromKey('InputsPath')
  refDataOPath  = INPUT.GetCharacterFromKey('OutputsPath')
  if 'SurrogateOnRefPoints' in keys:
    surOnRef = INPUT.GetLogicalFromKey('SurrogateOnRefPoints')
  else:
    surOnRef = False
  if "PlotKDE" in keys:
    refKDE = INPUT.GetLogicalFromKey('PlotKDE')
  else:
    refKDE = False
  if "SmoothenKDE" in keys:
    refKDESmooth = INPUT.GetLogicalFromKey('SmoothenKDE')
  else:
    refKDESmooth = False
  if "Preprocessing" in keys:
    isRefPreProc = True
    refPreProc = INPUT.GetCharacterFromKey('Preprocessing')
  else:
    isRefPreProc = False
except:
  pass
# ----- Loading Experiment -----
if isReference: #TODO: Check if a more relevant alternative would not be to load the reference experiment
  myExp = Experiment(empty=True, verbosity=1)
  myExp.Load(ExpeFile, prefix=ExpePath)

##############################
#      LOADING SURROGATE     #
##############################
if surType == 'Surrogate_PCE':
  mySur = PCE(empty=True, verbosity=1)
elif surType == 'Surrogate_Kriging':
  mySur = Kriging(empty=True, verbosity=1)
elif surType == 'Surrogate_1DSpline':
  mySur = OneDSpline(empty=True, verbosity=1) 
else:
  print('Surrogate Type "' + surType + '" not implemented.')
  quit()
mySur.Load(surFile, prefix=surPath)

##############################
#  RESAMPLING ORIGINAL DATA  #
##############################
if surOnOri:
  # ----- Loading data used to generate the surrogate -----
  datFile     = h5py.File(surFile, 'r')
  xOriEval    = datFile[surDataOPath][:]
  xOriSample  = datFile[surDataIPath][:,:]
  datFile.close()
  nOri        = len(xOriEval)
  xOriSurEval = np.zeros(nOri)

  # ----- Preprocessing and resampling original points -----
  if isSurPreProc:
    INPUT.GetBlock(surPreProc)
    xOriSample, dummy = Preprocessor(xOriSample, None, INPUT)
  if mySur.objectType == 'Surrogate_1DSpline':
    xOriSurEval = mySur.Evaluate(xOriSample)
  else:
    for i in range(nOri):
      xOriSurEval[i] = mySur.Evaluate(xOriSample[i])[0]
  if isSurPreProc:
    INPUT.GetBlock(surPreProc)
    dummy, xOriSurEval = Preprocessor(None, xOriSurEval, INPUT)
  if isOriPreProc:
    INPUT.GetBlock(oriPreProc)
    dummy, xOriEval = Preprocessor(None, xOriEval, INPUT)


##############################
# RESAMPLING VALIDATION DATA #
##############################
if isValidation:
  # ----- Loading validation data -----
  datFile     = h5py.File(valFile, 'r')
  xValEval    = datFile[valDataOPath][:]
  xValSample  = datFile[valDataIPath][:,:]
  datFile.close()
  nVal        = len(xValEval)
  xValSurEval = np.zeros(nVal)

  # ----- Preprocessing and resampling validation points -----
  if isSurPreProc:
    INPUT.GetBlock(surPreProc)
    xValSample, dummy = Preprocessor(xValSample, None, INPUT)
  if mySur.objectType == 'Surrogate_1DSpline':
    xValSurEval = mySur.Evaluate(xValSample)
  else:
    for i in range(nVal):
      xValSurEval[i] = mySur.Evaluate(xValSample[i])[0]
  if isSurPreProc:
    INPUT.GetBlock(surPreProc)
    dummy, xValSurEval = Preprocessor(None, xValSurEval, INPUT)
  if isValPreProc:
    INPUT.GetBlock(valPreProc)
    dummy, xValEval = Preprocessor(None, xValEval, INPUT)


##############################
# RESAMPLING REFERENCE DATA  #
##############################
if isReference:
  # ----- Loading reference data -----
  datFile     = h5py.File(refFile, 'r')
  xRefEval    = datFile[refDataOPath][:]
  xRefSample  = datFile[refDataIPath][:,:]
  datFile.close()
  nRef        = len(xRefEval)
  xRefSurEval = np.zeros(nRef)

  if surOnRef:
    # ----- Preprocessing and resampling reference points -----
    if isSurPreProc:
      INPUT.GetBlock(surPreProc)
      xRefSample, dummy = Preprocessor(xRefSample, None, INPUT)
    if mySur.objectType == 'Surrogate_1DSpline':
      xRefSurEval = mySur.Evaluate(xRefSample)
    else:
      for i in range(nRef):
        xRefSurEval[i] = mySur.Evaluate(xRefSample[i])[0]
    if isSurPreProc:
      INPUT.GetBlock(surPreProc)
      dummy, xRefSurEval = Preprocessor(None, xRefSurEval, INPUT)
    if isRefPreProc:
      INPUT.GetBlock(refPreProc)
      dummy, xRefEval = Preprocessor(None, xRefEval, INPUT)
  else:
    # ----- Preprocessing and sampling on new MC points -----
    distribution  = myExp.distribution
    nSample       = 100000
    xRefSurEval = np.zeros(nSample)
    MCDrawer = MC(None, distribution, nSample, firstSample=0, verbosity=1, fixSeed=False)
    MCDrawer.Draw()
    xNewSample = MCDrawer.xPoint
    if isSurPreProc:
      INPUT.GetBlock(surPreProc)
      xNewSample, dummy = Preprocessor(xNewSample, None, INPUT)
    #if mySur.objectType == 'Surrogate_1DSpline':
    #  xRefSurEval = mySur.Evaluate(xNewSample)
    #else:
    for i in range(nRef):
      xRefSurEval[i] = mySur.Evaluate(xNewSample[i])[0]
    if isSurPreProc:
      INPUT.GetBlock(surPreProc)
      dummy, xRefSurEval = Preprocessor(None, xRefSurEval, INPUT)


##############################
#       GLOBAL METRICS       #
##############################
if isValidation:
  # cf Khalil 2015
  NRMSE = norm2(  xValEval - xValSurEval) / norm2(  xValEval)
  NMAE  = norminf(xValEval - xValSurEval) / norminf(xValEval)
  prt('NRMSE =                  ' + str(round(NRMSE*100,3)) + ' %', 'None', True)
  prt('NMAE  =                  ' + str(round(NMAE*100,3)) + ' %', 'None', True)
#if isReference: #DEBUG
#  prt('Distance to reference =  ' + str(wdist(xRefEval, xRefSurEval)), "None", True)


############################## #DEBUG
#        SUMMARY PLOT        #
##############################
minX = xOriEval.min()
maxX = xOriEval.max()
if isValidation:
  minX = min(minX, xValEval.min())
  maxX = max(maxX, xValEval.max())
if isReference:
  minX = min(minX, xRefEval.min())
  maxX = max(maxX, xRefEval.max())
x1 = np.linspace(minX, maxX, 1000)

fig = plt.figure(1, figsize=(7,5))
ax1 = fig.add_subplot(111)
if isReference and surOnRef:
  print("plotting reference")
  ax1.scatter(np.log10(np.exp(xRefEval)), np.log10(np.exp(xRefSurEval)), color='k', alpha = 1.0, marker=',', s=1)#, label='Reference samples')
#if surOnOri:
#  ax1.scatter(xOriEval, xOriSurEval, color='r', alpha = 0.5, marker='o', s=1, label='Original samples')
if isValidation:
  ax1.scatter(xValEval, xValSurEval, color='b', alpha = 0.5, marker='x', s=1, label='Validation samples')
ax1.plot(np.log10(np.exp(x1)),np.log10(np.exp(x1)), color='#26c6da', linewidth=3, label='Ideal')
ax1.set_xlabel(r'Model evaluation(log$_{10}(\tau)$)')
ax1.set_ylabel(r'Surrogate evaluation(log$_{10}(\tau)$)')
#ax1.set_title("Summary plot")
plt.legend(loc='best')
#pp.CorrectAxSize(ax1)
plt.savefig(sumFile, bbox_inches='tight')


############################
#        MODEL PLOT        #
############################
if modFile != "":
  if mySur.dimension == 1:
    minX = xOriSample.min()
    maxX = xOriSample.max()
    if isValidation:
      minX = min(minX, xValSample.min())
      maxX = max(maxX, xValSample.max())
    if isReference and surOnRef:
      minX = min(minX, xRefSample.min())
      maxX = max(maxX, xRefSample.max())
    nPoint = 1000
    x2 = np.linspace(minX, maxX, nPoint)
    y2 = np.zeros_like(x2)
    if mySur.objectType == 'Surrogate_1DSpline':
      y2 = mySur.Evaluate(x2)
    else:
      for i in range(nPoint):
        y2[i] = mySur.Evaluate(x2[i])[0]
    if isSurPreProc:
      INPUT.GetBlock(surPreProc)
      dummy, y2 = Preprocessor(None, y2, INPUT)

    fig = plt.figure(1, figsize=(7, 5))
    ax1 = fig.add_subplot(111)
    ax1.plot(x2,y2, color='k', label=r'Surrogate')
    if isReference and surOnRef:
      ax1.scatter(xRefSample, xRefEval, color='gray', alpha=0.5, marker=',', s=1, label='Reference samples')
    if isValidation:
      ax1.scatter(xValSample, xValEval, color='b', alpha=0.5, marker='x', s=1, label='Validation samples')
    #ax1.scatter(xOriSample[:129], xOriEval[:129], color='r', alpha=0.5, marker='x', s=20, label='Construction samples')

    plt.legend(loc='best')
    #pp.CorrectAxSize(ax1)
    for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
      item.set_fontsize(20)
    plt.savefig(modFile, bbox_inches='tight')
    #plt.show()

  elif mySur.dimension == 2:

    from mpl_toolkits.mplot3d import axes3d
    from matplotlib import cm
    # DEBUG. TODO: better thingy here
    xMin = xRefSample[:,0].min()
    xMax = xRefSample[:,0].max()
    yMin = xRefSample[:,1].min()
    yMax = xRefSample[:,1].max()
    nPt = 100
    x = np.linspace(xMin, xMax, nPt)
    y = np.linspace(yMin, yMax, nPt)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    for i in range(nPt):
      for j in range(nPt):
        Z[i, j] = mySur.Evaluate(np.array([X[i, j], Y[i, j]]))[0]

    norm = plt.Normalize(Z.min(), Z.max())
    colors = cm.hot(norm(Z))
    rcount, ccount, _ = colors.shape

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(X, Y, Z, cmap=cm.coolwarm, alpha=0.3, label='Surrogate')
    #ax.scatter(xRefSample[:,0], xRefSample[:,1], xRefEval[:], marker=',', s=1, color='k', label='Reference samples')
    ax.scatter(xRefSample[:,0], xRefSample[:,1], np.log(xRefEval[:]), marker=',', s=1, color='k', label='Reference samples') #DEBUG
    #ax.scatter(xOriSample[:129,0], xOriSample[:129,1], np.log(xOriEval[:129]), marker='x', s=40, color='r', label='Original samples') #DEBUG
    #ax.view_init(azim=100, elev=9)
    ax.view_init(azim=0, elev=90)
    #ax.view_init(azim=-80, elev=9)
    ax.set_ylabel(r'$\hat{T}$', labelpad=15)
    ax.set_xlabel(r'$\xi_{16}$', labelpad=15)
    #ax.set_xlabel(r'$' + myExp.xVarID[0] + '$', labelpad=15)
    #ax.set_ylabel(r'$' + myExp.xVarID[1] + '$', labelpad=15)
    ax.set_zlabel(r'QoI', labelpad=15)
    #pp.CorrectAxSize(ax)
    plt.legend(loc='best')
    #plt.savefig(modFile, bbox_inches='tight') #DEBUG
    plt.show()


  else:
    # Nothing to do
    pass

# ---------------- END DEBUG --------------------


#plt.cla()
#plt.clf()
#plt.close()



###############################
##          PDF PLOT          #
###############################
#if isReference:
#  nColumn = 50 # TODO: Better handling of the num of bins
#  fig = plt.figure(figsize=(7,5))
#  ax = fig.add_subplot(111)
#  if refKDE:
#    refPDF = pp.KDE(xRefEval)
#    refPDF.OptimizeBandwidth(verySmooth=refKDESmooth)
#    xRef, yRef = refPDF.GetDataPlot()
#    ax.plot(xRef, yRef, color='k', label='PDF (reference)')
#  else:
#    n, bins, patches = ax.hist(xRefEval, nColumn, density=True, facecolor='dimgrey', alpha=0.33, label='PDF (reference)')
#  if surKDE:
#    surPDF = pp.KDE(xRefSurEval)
#    surPDF.OptimizeBandwidth(verySmooth=surKDESmooth)
#    xSur, ySur = surPDF.GetDataPlot()
#    ax.plot(xSur, ySur, color='b', linestyle=':', label='PDF (surrogate)')
#  else:
#    n, bins, patches = ax.hist(xRefSurEval, nColumn, density=True, facecolor='b', alpha=0.33, label='PDF (surrogate)')
#
#  #ax.set_xlabel('QoI')
#  #ax.set_ylabel('pdf(QoI)')
#  ax.set_xlabel(r'log($\tau$)')
#  ax.set_ylabel(r'pdf(log($\tau$))')
#  ax.set_xlim([xRefEval.min(), xRefEval.max()])
#  ax.set_title('PDF comparison')
#  ax.legend(loc='best')
#  pp.CorrectAxSize(ax)
#  plt.savefig(pdfFile, bbox_inches='tight')


#plt.clf()
#plt.cla()
#plt.close()
#------------------ END TEST

prt("THE END", 'green', True)
quit()
