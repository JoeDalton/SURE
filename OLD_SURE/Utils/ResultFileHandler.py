import h5py
import sys
import copy
import numpy as np
from I_O import prt

#TODO: Maybe this would need to be moved elsewhere


def createNewResultFile(resultFile, resultQoIPath, resultInputPath, nSample, nVar):
  # ----- Overwriting, result file -----
  result = h5py.File(resultFile, 'w')
  #TODO: Clean this with proper groups and datasets in the Header
  #result.create_group("Header")
  #result.create_dataset("Header/simuType", data="")
  #result.create_dataset("Header/QoI_Detection", data="")
  #result.create_dataset("Header/QoI_Value", data="")
  #result.create_dataset("Header/QoI_PhysQuantity", data="")
  result.create_group("Data")
  result.create_dataset(resultQoIPath, data=np.zeros(nSample))
  result.create_dataset(resultInputPath, data=np.zeros((nSample, nVar)))
  result.close()




def HandleResultFile(resultFile, resultQoIPath, resultInputPath, nSample, nVar, isContinued):
  #TODO: Better error information
  # ----- handling result files -----
  willExit = False
  try:
    result = h5py.File(resultFile, 'r+')
    if not isContinued:
      inp = raw_input(resultFile + " already exists. Overwrite ? (y/n) ")
      if inp in ['Y', 'y', 'Yes', 'yes']:
        result.close()
        os.system("rm " + resultFile)
        createNewResultFile(resultFile, resultQoIPath, resultInputPath, nSample, nVar)
        firstSample = 0
      else:
        result.close()
        prt("Exiting", 'red', True)
        willExit = True
    else:
      firstSample = len(result[resultQoIPath][:])
      tempQoI = result[resultQoIPath][:]
      tempNormalInputVars = result[resultInputPath][:,:]
      del result[resultQoIPath]
      del result[resultInputPath]
      tempQoI2 = np.zeros(nSample + firstSample)
      tempQoI2[:firstSample] = tempQoI[:]
      tempNormalInputVars2 = np.zeros((nSample+firstSample, nVar))
      tempNormalInputVars2[:firstSample,:] = tempNormalInputVars[:,:]
      result.create_dataset(resultQoIPath, data=tempQoI2)
      result.create_dataset(resultInputPath, data=tempNormalInputVars2)
      result.close()
      
  except:
    if isContinued:
      prt(resultFile + " does not exist. Exiting...", 'red', True)
      willExit = True
    else:
      createNewResultFile(resultFile, resultQoIPath, resultInputPath, nSample, nVar)
      firstSample = 0
  if willExit:
    sys.exit()
  return firstSample
