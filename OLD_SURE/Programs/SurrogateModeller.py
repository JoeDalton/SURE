#!/usr/bin/env python2.7

import os
import sys
import h5py
import Link
import argparse
import openturns                                            as ot
from copy                               import copy         as cp
from Utils.I_O                          import prt, ErrorMsg
from Utils.I_O                          import InputHandler as IH
from Experiment.Experiment              import *
from SurrogateModelling.PCE             import *
from SurrogateModelling.OneDSpline      import *
from SurrogateModelling.Kriging         import *
from PoolUsefulScripts.OK_model         import *
from DimensionReduction.ActiveDirection import *

fileName = 'Programs/SurrogateModeller.py'

# ----- Parsing arguments -----
parser = argparse.ArgumentParser(description="SurrogateModeller")
parser.add_argument('-i',	  required=False,	type=str,  help = "Relative path of input file")
args = parser.parse_args()

# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'SurrogateModeller.input'
  try:
    INPUT = IH(inputFile, modify=False)
  except:
    prt('Input File does not exist', 'None', True)
    prt('  --> Copy an example in ' + inputFile, 'None', True)
    os.system("cp " + Link.sureHome + '/InputFiles/' + inputFile + ' .')
    sys.exit()

# ----- Experiment-related reading-----
INPUT.GetBlock('Experiment')
myExperiment = Experiment(block=INPUT, simulation=None)

      


# ----- Data-related reading -----
INPUT.GetBlock('Data')
resultFile  = INPUT.GetCharacterFromKey('dataFile')
qoiPath     = INPUT.GetCharacterFromKey('QoIPath')
inputPath   = INPUT.GetCharacterFromKey('InputPath')
keys        = INPUT.GetKeysForBlock()
#if "Transformation" in keys:
#  transform = INPUT.GetCharacterFromKey('Transformation')
#else:
#  transform = 'None'
# ----- Extract data from result file -----
result      = h5py.File(resultFile, 'r')
# Get evaluations
if "nSample" in keys: # TODO: Clean this by sending it in a preprocessing section ?
  nSample   = INPUT.GetIntegerFromKey('nSample')
else:
  nSample   = len(result[qoiPath][:])
#if transform == 'log':
#  xEval     = np.log(result[qoiPath][:nSample])
#elif transform == 'None':
xEval     = result[qoiPath][:nSample]
#else:
#  prt('', 'red', True)
#  prt('SurrogateModeller', 'red', True)
#  prt('Error: The supported data transformation are "log" and "None".', 'red', True)
#  prt('Exiting', 'red', True)
#  sys.exit() 
# Get sample points
xSample   = result[inputPath][:nSample,:]

# Get weights
xWeight = []
if 'WeightPath' in keys: # TODO: remove this and rather check if PCE and Projection method, then ok for weights
  prt('', 'red', True)
  inp = raw_input("Warning: Weights have been provided. Are you sure you want to use them ? (y/n) ")
  if inp in ['Y', 'y', 'Yes', 'yes']:
    weightPath  = INPUT.GetCharacterFromKey('WeightPath')
    xWeight     = result[weightPath][:nSample]
  else:
    prt('Exiting', 'red', True)
    sys.exit() 
result.close()














# ----- Preprocessing -----
# TODO: Put this in a nice and sleek object ?
# Advantages : 
# 1. Hide this whole big stuff in a few elegant lines
# 2. Easy to save cleanly in the result file
# 3. Easily reusable for the other programs
# Drawbacks:
# 1. Long and complicated to write
# 2. I'm too lazy to do it right now
# Conclusion:
# I'll do it later. Maybe. Wait for it...
# https://www.youtube.com/watch?v=VBlFHuCzPgY

# Determine preprocessings to do
try:
  INPUT.GetBlock('Preprocessing')
  isPreproc   = True
  keys        = INPUT.GetKeysForBlock()
  # Pre-processing of sampling points
  if "xInPreProc" in keys:
    temp = INPUT.GetCharacterFromKey('xInPreProc')
    if isinstance(temp, list):
      xInPreProc  = temp
    else:
      xInPreProc  = [temp]
  else:
    xInPreProc    = []
  # Pre-Processing of sampled QoIs
  if "xOutPreProc" in keys:
    temp = INPUT.GetCharacterFromKey('xOutPreProc')
    if isinstance(temp, list):
      xOutPreProc = temp
    else:
      xOutPreProc = [temp]
  else:
    xOutPreProc   = []
except:
    xInPreProc    = []
    xOutPreProc   = []

# Preprocess the QoIs # Must be done first so that an eventual active substace will be correctly computed
tempOut = cp(xEval)
for proc in xOutPreProc:
  INPUT.GetSubBlock(proc)
  preType       = INPUT.GetCharacterFromKey('Type')
  if preType == 'LogTransform':
    tempOut = np.log(tempOut)
  else:
    ErrorMsg('The supported QoI preprocessing are: "LogTransform".', fileName) # TODO: Accept more, for example the scalings, and perhaps also other analyticals if necessary
xProcEval = tempOut

# Preprocess the sample points
tempIn = cp(xSample)
for proc in xInPreProc:
  INPUT.GetSubBlock(proc)
  preType       = INPUT.GetCharacterFromKey('Type')
  if preType == 'ActiveDirection':
    saveFile      = INPUT.GetCharacterFromKey('SaveFile')
    savePath      = INPUT.GetCharacterFromKey('SavePath')
    myAD = ActiveDirection(block=INPUT, experiment=myExperiment) # Warning: not rigorous if it is not the first applied...
    myAD.Build(xSample, xEval)
    myAD.Save()
    myAD.Dump(fileName=saveFile, path=savePath)
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








# ----- Surrogate-related reading -----
INPUT.GetBlock('Surrogate_Model')
surrogateType = INPUT.GetCharacterFromKey('Type')
saveFile      = INPUT.GetCharacterFromKey('SaveFile')
savePath      = INPUT.GetCharacterFromKey('SavePath')

if surrogateType == 'PCE':
  mySurrogate = PCE(block=INPUT, experiment=myExperiment)
  if len(xWeight) == 0:
    mySurrogate.Build(xProcSample, xProcEval)
  else:
    mySurrogate.Build(xProcSample, xProcEval, xWeight=xWeight)
elif surrogateType == 'OneDSpline':
    surrogateID = INPUT.GetCharacterFromKey('ID')
    #mySurrogate = OneDSpline(ID=surrogateID, verbosity=1)
    mySurrogate = OneDSpline(block=INPUT, experiment=myExperiment)
    mySurrogate.Build(samples=xProcSample, evaluations=xProcEval)
elif surrogateType == 'Kriging':
    #mySurrogate = OK_model(block=INPUT, experiment=myExperiment)
    mySurrogate = Kriging(block=INPUT, experiment=myExperiment)
    mySurrogate.Build(xProcSample, xProcEval)
else:
  ErrorMsg('The supported surrogate types are "PCE", "Kriging" and "OneDSpline".', fileName)

# ----- Saving the surrogate -----
mySurrogate.Save()
mySurrogate.Dump(fileName=saveFile, path=savePath)





prt("THE END", 'green', True)
quit()
