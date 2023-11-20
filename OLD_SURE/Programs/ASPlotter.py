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
parser      = argparse.ArgumentParser(description="AS plotter")
parser.add_argument('-i',	          required=False,	type=str,   help = "Relative path of input file")
args        = parser.parse_args()

# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'ASPlotter.input'
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
INPUT.GetBlock('AS_Coefficients')
myFile      = INPUT.GetCharacterFromKey('LoadFile')
myPath      = INPUT.GetCharacterFromKey('LoadPath')
threshold   = INPUT.GetRealFromKey('PlotThreshold')
asFile      = h5py.File(myFile, 'r')
indices     = asFile[myPath][:]
asFile.close()

prt("Sensitivity indices loaded", 'green', True) #DEBUG
dimension   = len(indices)




#------------------ TEST: To move to Sensitivity plotter

# ----- Build a dictionary with all the indices
asDict = {}
# Add first-order indices
for i in range(dimension):
    myKey           = '$' + xVarName[i] + '$'
    myValue         = indices[i]
    asDict[myKey]   = myValue

# ----- Sort in decreasing order
xKey    = np.array(asDict.keys())
xVal    = np.array(asDict.values())
xValAbs = np.abs(np.array(asDict.values()))
mySort  = np.argsort(xValAbs)[::-1]
sKey    = xKey[mySort]
sVal    = xVal[mySort]

# ----- Keep only the ones summing to the chosen fraction of variance
mySum = 0.0
majI  = len(sKey)
for i in range(len(sKey)):
  if mySum < (1.0-threshold):
    mySum += sVal[i] * sVal[i]
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
colors=['r' if sVal[i]>0 else 'b' for i in range(majI)]
#ax1.scatter(range(majI), sVal[:majI] * sVal[:majI], color='k', marker='o', facecolors='none')
ax1.scatter(range(majI), sVal[:majI] * sVal[:majI], color=colors, marker='o', facecolors='none')
ax1.scatter([],[], color='b', marker='o', facecolors='none', label='Negative contribution')
ax1.scatter([],[], color='r', marker='o', facecolors='none', label='Positive contribution')

ax1.set_ylim([1e-2, 1e0])
ax1.set_xticks(range(majI))
#ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18) #DEBUG
#ax1.set_xticklabels([r'$T$', r'$A_{16}$', r'$A_{17}$', r'$A_{20}$', r'$A_{24}$', r'$A_{15}$', r'$A_{22}$', r'$A_{26}$'], rotation='vertical', fontsize=18)
ax1.set_xticklabels([r'$T$', r'$A_{16}$'], rotation='vertical', fontsize=18)

ax1.set_ylabel(r'Sensitivity indices (${w_{1,j}}^2$)')
plt.yscale('log')
plt.legend()
for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(18)

#plt.show()
plt.savefig('asplot.pdf', bbox_inches='tight') #DEBUG


prt("THE END", 'green', True)
quit()
