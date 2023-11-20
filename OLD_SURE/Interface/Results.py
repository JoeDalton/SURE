import h5py
import sys
import os
import math
from Utils.I_O import prt, ErrorMsg, WarningMsg
import numpy as np
from scipy.interpolate                          import RectBivariateSpline
from Numerics.Optimisation.Genetics.Genetics    import *
from Numerics.Optimisation.Genetics.Population  import *
from Numerics.Optimisation.Genetics.Eval        import MultiEval

fileName = 'Interface/Results.py'

def AGATH_0D_IDT(directory):
  #TODO: Add a test of valid h5 file for robustness
  #TODO: Maybe add an option to change the result file name ? Checking the input file for example ?
  solutionFile = directory + "/Solution.h5"
  try:
    sol = h5py.File(solutionFile, 'r')
    isValidFile = False
    if "Header" in sol:
      if "Main_IgnitionDelay" in sol["Header"]:
        isValidFile = True
  
    if not isValidFile:
      ErrorMsg('The file ' + solutionFile + ' is not a valid OD reactor solution.', fileName)

    QoI = float(sol["Header"]["Main_IgnitionDelay"][0])
    sol.close()
    return QoI
  except:
    WarningMsg(directory + " seems to have not converged !",fileName, True)
    return 0.0

def Simulation_To_Keep(directory):
  # Always return 0.0 so that the simulation folder will not be erased.
  # TODO: Change this when the "non convergence" flag is also changer
  return 0.0

def AGATH_1D_IDT(directory):
  # QoI is the time stamp of the last solution
  
  # List all files in the solution directory
  directory += "/Sols" 
  from os import listdir
  from os.path import isfile, join
  allFiles = [f for f in listdir(directory) if isfile(join(directory, f))]
  
# Auto-ignited solution file finishes by "_last.h5"
  interestFile = ''
  for myFile in allFiles:
    if "Solution_" in myFile and "_last.h5" in myFile:
      interestFile = directory + '/' + myFile
  try:
    sol = h5py.File(interestFile, 'r')
    isValidFile = False
    if "Header" in sol:
      if "Time" in sol["Header"]:
        isValidFile = True
    if not isValidFile:
      ErrorMsg('The file ' + interestFile + ' is not a valid 1D reactor solution.', fileName)

    QoI = float(sol["Header"]["Time"][0])
    sol.close()
    return QoI
  except:
    # No corresponding file : Sample has not auto-ignited 
    return -1.0
