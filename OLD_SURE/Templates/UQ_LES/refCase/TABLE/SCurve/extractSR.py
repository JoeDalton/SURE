#!/usr/bin/env python2.7
from pylab import *
import h5py
import sys
import argparse
import string
import os
import numpy as np
import Link
from Interface.Modifications import ModifyInputFile


Files    = sort(os.listdir("./Sols_Scurve"))

#-----------------------------------------------------
# Read Solution Files
#-----------------------------------------------------
xSR = []
xMHR = []
for SolutFile in Files:
  f    = h5py.File("./Sols_Scurve/" + SolutFile, 'r')
  SR = float(f["Header/Strain_Rate"][0])
  MHR = max(f["Data/HeatRelease/Array"][:])
  xSR.append(SR)
  xMHR.append(MHR)
  f.close()
xSR = np.array(xSR)
xMHR = np.array(xMHR)

idxMax = np.argmax(xMHR)
SrCrit = np.min(xSR[idxMax:])


filePath = "../buildTable/RunExperiment.input"
dictPath = "Experiment/Variable_Strain_Rate/transformationValues"
toWrite  = "5 " + str(SrCrit)
ModifyInputFile(filePath, dictPath, toWrite)

