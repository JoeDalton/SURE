import h5py
import time
import numpy as np
import sys
import os
import subprocess
from pylab import *
import h5py


workDir = "./TABLE"


t0 = time.time()


# S-curve computation
with open(workDir + '/SCurve/outCre', 'w') as o:
  cmdList = ["CreateAgathSolution.e"]
  p = subprocess.Popen(cmdList, cwd=workDir + '/SCurve', stdout=o)
  p.wait()

with open(workDir + '/SCurve/outUns', 'w') as o:
  cmdList = ["UnsteadyDiffZFlame.e"]
  p = subprocess.Popen(cmdList, cwd=workDir + '/SCurve', stdout=o)
  p.wait()

with open(workDir + '/SCurve/outSte', 'w') as o:
  cmdList = ["DiffZFlame.e"]
  p = subprocess.Popen(cmdList, cwd=workDir + '/SCurve', stdout=o)
  p.wait()
print("First flamelet converged, initiating S-curve computation...")

with open(workDir + '/SCurve/outMas', 'w') as o:
  cmdList = ["Master_SteadyFlamelet.e"]
  p = subprocess.Popen(cmdList, cwd=workDir + '/SCurve', stdout=o)
  p.wait()
t1 = time.time()

p = subprocess.Popen("PlotParams.py --savefig SCurve.pdf --sol Sols_Scurve/* --var Strain_Rate MHeatRelease", cwd=workDir + '/SCurve', shell=True)
p.wait()
print("SCurve plotted, " + str(t1-t0)+"s")


# Second critical point extraction
cmdList = ["python","extractSR.py"]
p = subprocess.Popen(cmdList, cwd=workDir + '/SCurve')
p.wait()
t2 = time.time()
print("Second critical point extracted and communicated, " + str(t2-t1)+"s")


# Flamelets generation
p = subprocess.Popen('RunExperiment.py', cwd=workDir + '/buildTable')
p.wait()
t3 = time.time()
print("Flamelets computed, " + str(t3-t2)+"s")


# Table dressing
with open(workDir + '/SCurve/outCre', 'w') as o:
  cmdList = ["bash", "MasterTable.e"]
  p = subprocess.Popen(cmdList, cwd=workDir + '/buildTable', stdout=o)
  p.wait()
t4 = time.time()
print("Table dressed, " + str(t4-t3)+"s")


# Table conversion
with open(workDir + '/SCurve/outCre', 'w') as o:
  cmdList = ["bash", "ConvertTable.e"]
  p = subprocess.Popen(cmdList, cwd=workDir + '/buildTable', stdout=o)
  p.wait()
t5 = time.time()
print("Table converted, " + str(t5-t4)+"s")


# Table completion
Files = sort(os.listdir(workDir + "/SCurve/Sols_Scurve"))
xSR = []
xMHR = []
for SolutFile in Files:
  f    = h5py.File(workDir + "/SCurve/Sols_Scurve/" + SolutFile, 'r')
  SR = float(f["Header/Strain_Rate"][0])
  MHR = max(f["Data/HeatRelease/Array"][:])
  xSR.append(SR)
  xMHR.append(MHR)
  f.close()
xSR = np.array(xSR)
xMHR = np.array(xMHR)
idxMax = np.argmax(xMHR)
SrCrit = np.min(xSR[idxMax:])

cmdList = ["python", "FinishUFPVV2.py", "-t", "UFPV_AVBP.h5", "--dummySR", str(SrCrit)]
p = subprocess.Popen(cmdList, cwd=workDir + '/buildTable')
p.wait()
cmdList = ["mv", "UFPV_AVBP.h5", "../."]
p = subprocess.Popen(cmdList, cwd=workDir + '/buildTable')
p.wait()
t6 = time.time()
print("Table completed, " + str(t6-t5)+"s")


print("The End !")
