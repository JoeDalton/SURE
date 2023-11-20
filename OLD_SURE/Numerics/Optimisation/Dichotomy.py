import numpy as np
import math
import sys
from Utils.I_O import ErrorMsg

fileName = 'Numerics/Optimisation/Dichotomy.py'

################################################################
# Optimisation with a dichotomy algorithm. 
# Options are "linear" for "usual" functions and "geometric" for functions with exponential behavior
############################################################
def DichotomyOptim(func, lb, ub, nInit, rtol, dichoType='linear'):
  # ---------------------------------------------------------------------
  # Inner function
  # ---------------------------------------------------------------------
  def Dichotomise(func, inSamples, inEvals, dichoType):
    # Find the most interesting sample
    mostInt       = np.argmin(inEvals)
    # Create new samples around the currently most interesting
    if mostInt == 0 or mostInt == len(inSamples)-1:
      newSamples    = np.zeros(1)
      newEvals      = np.zeros(1)
    else:
      newSamples    = np.zeros(2)
      newEvals      = np.zeros(2)
    if dichoType == 'geometric':
      if mostInt == 0:
        newSamples[0] = math.sqrt(inSamples[mostInt] * inSamples[mostInt+1])
      elif mostInt == len(inSamples)-1:
        newSamples[0] = math.sqrt(inSamples[mostInt] * inSamples[mostInt-1])
      else:
        newSamples[0] = math.sqrt(inSamples[mostInt] * inSamples[mostInt+1])
        newSamples[1] = math.sqrt(inSamples[mostInt] * inSamples[mostInt-1])
    elif dichoType == 'linear':
      if mostInt == 0:
        newSamples[0] = 0.5 * (inSamples[mostInt] + inSamples[mostInt+1])
      elif mostInt == len(inSamples)-1:
        newSamples[0] = 0.5 * (inSamples[mostInt] + inSamples[mostInt-1])
      else:
        newSamples[0] = 0.5 * (inSamples[mostInt] + inSamples[mostInt+1])
        newSamples[1] = 0.5 * (inSamples[mostInt] + inSamples[mostInt-1])
    else:
      ErrorMsg('Error: Please provide a valid dichoType argument. Valid arguments are "geometric" and "linear".', fileName)
    if mostInt == 0 or mostInt == len(inSamples)-1:
      newEvals[0]   = func(newSamples[0])
    else:
      newEvals[0]   = func(newSamples[0])
      newEvals[1]   = func(newSamples[1])
    # Add the new samples to the existing ones
    outSamples    = np.append(inSamples, newSamples)
    outEvals      = np.append(inEvals, newEvals)
    # Sort the lists in increasing smoothing factor (samples list)
    outSort       = outSamples.argsort()
    outSamples    = outSamples[outSort]
    outEvals      = outEvals[outSort]
    return outSamples, outEvals
  # ---------------------------------------------------------------------
  
  # Phase 1.0: Exploration on nInit samples spaced evenly on a log scale
  mySamples           = np.geomspace(lb, ub, nInit) 
  myEvals             = np.zeros(nInit)
  for it in range(nInit):
    myEvals[it]       = func(mySamples[it])
  # Phase 1.1: Find the most interesting value to initialize the dichotomy
  mySort              = myEvals.argsort()
  # Phase 2 : Dichotomy algorithm
  while True:
    mySamples, myEvals = Dichotomise(func, mySamples, myEvals, dichoType)
    mySort             = myEvals.argsort()
    test               = myEvals[mySort]
    residual           = (test[1]-test[0])/test[1]
    if residual<rtol:
      break

  mySamples   = mySamples[mySort]
  myBestInput = mySamples[0]
  return myBestInput
