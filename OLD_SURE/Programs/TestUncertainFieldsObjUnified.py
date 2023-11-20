#!/usr/bin/env python2.7

from __future__ import print_function
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc
from matplotlib.colors import ListedColormap
import h5py
import Link
from DimensionReduction.Preprocessing import *
from DimensionReduction.PCA import *
from DimensionReduction.ActiveDirection import *
from SurrogateModelling.Kriging import *
from SurrogateModelling.PCE import *
from SurrogateModelling.OneDSpline import *
from Experiment.Experiment import *
from Sampling.QMC import *
from Utils.I_O import prt, LoadCSV
import Utils.PrettyPlots as pp
import matplotlib as mpl #TODO: move this in a special module (Utils.I_O ?). It works and would purify the programs!
font = {'family' : 'sans',
   'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)

from PoolUsefulScripts.OK_model import *
from Numerics.Statistics import *
from Experiment.VarTransformation import *
from Numerics.Algebra import NMAE, NRMSE
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import os

################################################################
# STEP 0: Important parameters
############################################################
nSample       = 20       # Number of observations, or samples, considered in the study
nFeature = 78750
threshold = 2.3e-5 #
D =  0.00457
xstep = 10
rstep = 2

epsilon = 0.01 # Fraction of variance left on the side of the road


#
variable = 'OH'
#variable = 'z_out'
#variable = 'T'
#variable = 'H2'
cLab = r'$' + variable + '$ [-]'


reconstructedSample = 0



xVar = ['OH', 'O2', 'H2', 'H2O', 'T', 'z_out']
#xVar = ['OH', 'O2', 'H2', 'H2O', 'T']
nField = len(xVar)


################################################################
# STEP 1: Loading all the samples
############################################################
print("")
print("--------------------------------------------------------------")
prt('Loading samples', 'yellow', True)

solFile = '../results/UQ_LES_26.h5'
sol = h5py.File(solFile, 'r')
for j in range(nField):
  var = xVar[j]
  tempxVect = np.zeros((nFeature, nSample))
  for i in range(nSample):
    solPath     = 'Fields/sample_' + str(i) + '/' + var + '/'
    resr        = sol[solPath + 'resr'][()]
    resx        = sol[solPath + 'resx'][()]
    rlim        = sol[solPath + 'rlim'][:]
    xlim        = sol[solPath + 'xlim'][:]
    xEval       = sol[solPath + 'Array'][:]
    tempxVect[:,i]  = xEval[:]
  xSample = sol["Data/Inputs"][:nSample,:]
  xTransformParamTemp = sol['Experiment/Variables/T/xTransformationParam'][:]
  if j ==0:
    xVect = tempxVect.copy()
  else:
    xVect = np.concatenate([xVect, tempxVect], axis=0)
sol.close()










print("")
print("--------------------------------------------------------------")
prt('Building field emulator...', 'yellow', True)

from FieldEmulator import MultipleFieldEmulator

myEmulator=MultipleFieldEmulator(ID=variable, nField=nField)
myEmulator.Build(xSample, xVect)


print("")
print("--------------------------------------------------------------")
prt('Resampling experiment...', 'yellow', True)

#############################################
#   Reconstructing the validation samples   #
#############################################
nReSample = 11
xQ        = np.zeros((nReSample, nFeature))
ExpeFile  = solFile
ExpePath  = 'Experiment'
myExp = Experiment(empty=True, verbosity=1)
myExp.Load(ExpeFile, prefix=ExpePath)
distribution  = myExp.distribution
MCDrawer      = QMC(None, distribution, nReSample, firstSample=0, verbosity=1, sequence='Sobol')
MCDrawer.Draw()
xReSample       = MCDrawer.xPoint

xQ = myEmulator.Evaluate_Group(xReSample)


print("")
print("--------------------------------------------------------------")
prt('Evaluating error...', 'yellow', True)
from CabraError import *
xErr = np.zeros(nReSample)
for i in range(nReSample):
  err = ComputeCabraError(xQ[i,:], resr, resx, rlim, xlim)
  print (err)
  xErr[i] = err


print("")
print("--------------------------------------------------------------")
prt('Building surrogate...', 'yellow', True)
myPCE = PCE(ID='PCE', verbosity=1, distribution=ot.ComposedDistribution([ot.Normal(0.0,1.0), ot.Uniform(0.0,1.0)]))
myPCE.verbosity = True
myPCE.DefinePolyBasis(quasiNorm=0.7, weights=None) 
myPCE.DefineTruncatureStrat(strategy="Fixed", maxTotDegree=4)
myPCE.DefineEvalStrat(evalStrategy="Regression", validStrategy="KFold")
myPCE.Build(xReSample, xErr)


print("")
print("--------------------------------------------------------------")
prt('Optimizing...', 'yellow', True)
#idxMin = np.argmax(xErr)
idxMin = 0
bestGuess = xReSample[idxMin,:]# sample with the lowest error

def fun(vect):
  return myPCE.Evaluate(vect)[0]

minimized = minimize(fun, bestGuess, bounds=[(-1.5,1.5), (0.0, 1.0)]).x

print(minimized)
print(myPCE.Evaluate(minimized))
  
  
  
  
quit()



































print("")
print("--------------------------------------------------------------")
prt('Building mean, upper and lower bound fields...', 'yellow', True)

confidenceInterval = 90

loP = (50.0 - confidenceInterval/2.0)
hiP = (50.0 + confidenceInterval/2.0)



meQ = np.clip(np.mean(xQ, axis=0), a_min=0, a_max=None)
noQ = np.clip(xVect[:,0], a_min=0, a_max=None)
loQ = np.clip(np.percentile(xQ, loP, axis=0), a_min=0, a_max=None)
hiQ = np.clip(np.percentile(xQ, hiP, axis=0), a_min=0, a_max=None)


print("")
print("--------------------------------------------------------------")
prt('Plotting...', 'yellow', True)


xq = np.linspace(xlim[0],xlim[1],resx)
rq = np.linspace(rlim[0],rlim[1],resr)


csvDir = '../../LETSFIXIT/ExpProfiles/Cabra2004Stripped'



xNX   = [1,9,11,14,26]


# Create fig
yLab = r'$' + variable + '$ [-]'
fig = plt.figure(figsize=(8.0, 7.0 * ( len(xNX) + 1 )) )

nPlot = len(xNX) + 1
xAx = []
for i in range (nPlot):
  xAx.append(plt.subplot(nPlot, 1, i+1))

#Print the crap




ax = xAx[0]

meF = np.reshape(meQ, (resr, resx))
noF = np.reshape(noQ, (resr, resx))
loF = np.reshape(loQ, (resr, resx))
hiF = np.reshape(hiQ, (resr, resx))

meF[np.isnan(meF)] = 0
loF[np.isnan(loF)] = 0
hiF[np.isnan(hiF)] = 0

ax.plot(xq, meF[0,:], color='k', label='Expectation')
ax.plot(xq, noF[0,:], color='r', label='Nominal')
ax.fill_between(xq, hiF[0,:], loF[0,:], color='k', alpha=0.5, label=r''+ str(confidenceInterval) + '\% Confidence interval')
#ax.plot(xq, loF[0,:], color='k', linestyle='-')
#ax.plot(xq, hiF[0,:], color='k', linestyle='-')
ax.legend(fontsize=12)



for i in range(1,len(xNX)+1):
  normAxialDistance = xNX[i-1]
  ax = xAx[i]
  meP = []
  noP = []
  loP = []
  hiP = []
  for j in range(resr):
    meI = interp1d(xq, meF[j,:], kind='cubic')
    meP.append(meI(normAxialDistance*D))

    noI = interp1d(xq, noF[j,:], kind='cubic')
    noP.append(noI(normAxialDistance*D))

    loI = interp1d(xq, loF[j,:], kind='cubic')
    loP.append(loI(normAxialDistance*D))

    hiI = interp1d(xq, hiF[j,:], kind='cubic')
    hiP.append(hiI(normAxialDistance*D))


  ax.plot(rq, meP, color='k', label='Expectation')
  ax.plot(rq, noP, color='r', label='Nominal')
  ax.fill_between(rq, hiP, loP, color='k', alpha=0.5, label=r''+ str(confidenceInterval) + '\% Confidence interval')
  ax.legend(fontsize=12)







# Decorate the crap
for i in range(len(xAx)):
  if i == 0:
    maxX = 50.0*D
    xticks = np.arange(0.0, maxX, rstep * D * 2.5)
    xAx[i].set_xlabel('x/D')
    xAx[i].set_title('Axial')
  else:
    maxX = 10.0*D
    xticks = np.arange(0.0, maxX, rstep * D)
    xAx[i].set_xlabel('r/D')
    xAx[i].set_title('X = ' + str(xNX[i-1]) + ' D')
  xtickslabels = np.round(xticks/D)
  xAx[i].set_aspect('auto')
  xAx[i].set_xlim([0.0, maxX])
  xAx[i].set_xticks(xticks)
  xAx[i].set_xticklabels(xtickslabels)
  xAx[i].set_ylabel(yLab)
  pp.CorrectAxSize(xAx[i])
pp.AdjustSubplotSpacings(doubleVert=True)


# Add exp data from csvs
if variable=="z_out":
  csvPrefx = csvDir + '/' + 'Z' + '_'
else:
  csvPrefx = csvDir + '/' + variable + '_'
csvSufx = '.csv'

# Axial
csvName = csvPrefx + 'Axial' + csvSufx
x,y = LoadCSV(csvName)
xAx[0].scatter(x*D,y, label='Exp', edgecolors='k',facecolors='none', marker='o')
xAx[0].legend(fontsize=12)

# Radial
for i in range(len(xNX)):
  csvName = csvPrefx + 'z:d=' + str(xNX[i]) + csvSufx
  x,y = LoadCSV(csvName)
  xAx[i+1].scatter(np.abs(x/1000.0),y, label='Exp', edgecolors='k',facecolors='none', marker='o')
  xAx[i+1].legend(fontsize=12)

#plt.savefig('MasterUncertainProfile_X_' + variable + '.pdf', bbox_inches='tight', pad_inches=0)
plt.show()


os.system('say -v "Daniel" "It is done. This message will self-destruct in 5. 4. 3. 2. 1. Boum"')
