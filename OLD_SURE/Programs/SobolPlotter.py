#!/usr/bin/env python2.7

import os
import sys
import h5py
import Link
import argparse
from Utils.I_O                import prt
from Utils.I_O                import InputHandler as IH
from Experiment.Experiment import Experiment as XP
import numpy as np

import matplotlib.pyplot as plt

# ----- Parsing arguments -----
parser      = argparse.ArgumentParser(description="SobolIndicesFinder")
parser.add_argument('-i',	          required=False,	type=str,   help = "Relative path of input file")
args        = parser.parse_args()

# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'SobolPlotter.input'
  try:
    INPUT = IH(inputFile, modify=False)
  except:
    prt('Input File does not exist', 'None', True)
    prt('  --> Copy an example in ' + inputFile, 'None', True)
    os.system("cp " + Link.sureHome + '/InputFiles/' + inputFile + ' .')
    sys.exit()

# ----- Experiment-related reading-----
INPUT.GetBlock('Experiment')
myExperiment  = XP(block=INPUT, simulation=None)
xVarName      = myExperiment.xVarID

# ----- Sobol-related reading-----
INPUT.GetBlock('Sobol_Indices')
myFile      = INPUT.GetCharacterFromKey('LoadFile')
myPath      = INPUT.GetCharacterFromKey('LoadPath')
threshold   = INPUT.GetRealFromKey('PlotThreshold')
firstPath   = myPath + '/first_order'
secondPath  = myPath + '/second_order'
totalPath   = myPath + '/total_order'
sobFile     = h5py.File(myFile, 'r')
firstOrder  = sobFile[firstPath][:]
secondOrder = sobFile[secondPath][:]
totalOrder  = sobFile[totalPath][:]
sobFile.close()

prt("Sensitivity indices loaded", 'green', True) #DEBUG
dimension   = len(firstOrder)




#------------------ TEST: To move to Sensitivity plotter

# ----- Build a dictionary with all the indices
sobolDict = {}
# Add first-order indices
for i in range(dimension):
    myKey           = '$' + xVarName[i] + '$'
    myValue         = firstOrder[i]
    sobolDict[myKey]  = myValue
# Add second-order indices
for i in range(dimension):
  for j in range(i):
    myKey           = '$' + xVarName[i] + ' \& ' + xVarName[j] + '$'
    myValue         = secondOrder[i,j]
    sobolDict[myKey]  = myValue

# ----- Sort in decreasing order
xKey    = np.array(sobolDict.keys())
xVal    = np.array(sobolDict.values())
mySort  = np.argsort(xVal)[::-1]
sKey    = xKey[mySort]
sVal    = xVal[mySort]

# ----- Keep only the ones summing to the chosen fraction of variance
mySum = 0.0
majI  = len(sKey)
for i in range(len(sKey)):
  if mySum < (1.0-threshold):
    mySum += sVal[i]
  else:
    majI = i
    break

# ----- Plot the significant indices in decreasing order

import matplotlib as mpl
font = {'family' : 'sans',
    'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)

#Vanilla plot
fig, ax1 = plt.subplots(1,1)
ax1.scatter(range(majI), sVal[:majI], color='k', marker='o', facecolors='none')

ax1.set_ylim([1e-4, 1e0])
ax1.set_xticks(range(majI))
ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

ax1.set_ylabel('Sobol index')
plt.yscale('log')

for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(20)

#plt.show()
plt.savefig('sobolplot.pdf', bbox_inches='tight') #DEBUG


prt("THE END", 'green', True)
quit()
