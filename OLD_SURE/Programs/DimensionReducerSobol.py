#!/usr/bin/env python2.7

import os
import sys
import h5py
import Link
import argparse
import openturns                                  as ot
from Utils.I_O                import prt
from Utils.I_O                import InputHandler as IH
from Utils.I_O                import DumpObject
from SurrogateModelling.PCE   import *
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

# Experiment-related
INPUT.GetBlock('Parameters')
inFile    = INPUT.GetCharacterFromKey('PCEFile')
inPath    = INPUT.GetCharacterFromKey('PCEPath')
if shouldWrite:
  outFile   = INPUT.GetCharacterFromKey('OutputFile')
  outPath   = INPUT.GetCharacterFromKey('OutputPath')


myPCE = PCE(empty=True, verbosity=1) 
myPCE.Load(inFile, prefix=inPath)

# Retrieve var names
PCEFile = h5py.File(inFile, 'r')
temp = PCEFile['Experiment/Variable_IDs'][:]
PCEFile.close()
xVarName = []
for var in temp:
  xVarName.append(r'$' + var + '$')


sensAnalysis  = ot.FunctionalChaosSobolIndices(myPCE.result)

firstOrder    = [sensAnalysis.getSobolIndex(i) for i in range(myPCE.dimension)]
totalOrder    = [sensAnalysis.getSobolTotalIndex(i) for i in range(myPCE.dimension)]

secondOrder   = np.zeros((myPCE.dimension, myPCE.dimension))
for i in range(myPCE.dimension):
  for j in range(i):
    secondOrder[i,j] = sensAnalysis.getSobolIndex([i, j])

if shouldWrite:
  toDump = {}
  toDump['first_order']   = firstOrder
  toDump['total_order']   = totalOrder
  toDump['second_order']  = secondOrder
  #DumpObject(toDump, outFile, outPath) #DEBUG


print("Sensitivity indices found") #DEBUG



#------------------ TEST: To move to Sensitivity plotter

# ----- Build a dictionary with all the indices
sobolDict = {}
# Add first-order indices
for i in range(myPCE.dimension):
    myKey           = xVarName[i]
    myValue         = firstOrder[i]
    sobolDict[myKey]  = myValue
# Add second-order indices
for i in range(myPCE.dimension):
  for j in range(i):
    myKey           = xVarName[i] + ' \& ' + xVarName[j]
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
  if mySum < (1.0-ignoredVar):
    mySum += sVal[i]
  else:
    majI = i
    break

# ----- Plot the significant indices in decresaing order

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
#ax1.errorbar(range(34), firstOrder, yerr=firstError, fmt='x', ecolor='k', ls='none') #, marker='+', label=r'first order indices')
#ax1.errorbar(range(34), totalOrder, yerr=totalError, ecolor='r') #, marker='o', facecolors='none', label=r'total order indices')

ax1.set_ylim([1e-4, 1e0])
ax1.set_xticks(range(majI))
ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

ax1.set_ylabel(r'Sobol index')
plt.yscale('log')

for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(18)

#plt.show()
plt.savefig('sobolplot.pdf', bbox_inches='tight') #DEBUG



#TEST2
from matplotlib.patches import Rectangle
#Keeping only first var
fig, ax1 = plt.subplots(1,1)
ax1.scatter(range(majI), sVal[:majI], color='k', marker='o', facecolors='none')
#ax1.errorbar(range(34), firstOrder, yerr=firstError, fmt='x', ecolor='k', ls='none') #, marker='+', label=r'first order indices')
#ax1.errorbar(range(34), totalOrder, yerr=totalError, ecolor='r') #, marker='o', facecolors='none', label=r'total order indices')

ax1.set_ylim([1e-4, 1e0])
ax1.set_xticks(range(majI))
ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

ax1.set_ylabel(r'Sobol index')
plt.yscale('log')

for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(18)

rect = Rectangle((0.5,1e-4),19,1e4,linewidth=1,edgecolor='none',facecolor='grey', alpha = 0.7)
ax1.add_patch(rect)
plt.vlines(0.5, 1e-4, 1, colors='r')

#plt.show()
plt.savefig('Sobol1.pdf', bbox_inches='tight') #DEBUG

#Keeping only 2 vars
fig, ax1 = plt.subplots(1,1)
ax1.scatter(range(majI), sVal[:majI], color='k', marker='o', facecolors='none')
#ax1.errorbar(range(34), firstOrder, yerr=firstError, fmt='x', ecolor='k', ls='none') #, marker='+', label=r'first order indices')
#ax1.errorbar(range(34), totalOrder, yerr=totalError, ecolor='r') #, marker='o', facecolors='none', label=r'total order indices')

ax1.set_ylim([1e-4, 1e0])
ax1.set_xticks(range(majI))
ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

ax1.set_ylabel(r'Sobol index')
plt.yscale('log')

for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(18)

rect = Rectangle((1.5,1e-4),18,1e4,linewidth=1,edgecolor='none',facecolor='grey', alpha = 0.7)
ax1.add_patch(rect)
plt.vlines(1.5, 1e-4, 1, colors='r')
#plt.show()
plt.savefig('Sobol2.pdf', bbox_inches='tight') #DEBUG


#Keeping only 3 vars
fig, ax1 = plt.subplots(1,1)
ax1.scatter(range(majI), sVal[:majI], color='k', marker='o', facecolors='none')
#ax1.errorbar(range(34), firstOrder, yerr=firstError, fmt='x', ecolor='k', ls='none') #, marker='+', label=r'first order indices')
#ax1.errorbar(range(34), totalOrder, yerr=totalError, ecolor='r') #, marker='o', facecolors='none', label=r'total order indices')

ax1.set_ylim([1e-4, 1e0])
ax1.set_xticks(range(majI))
ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

ax1.set_ylabel(r'Sobol index')
plt.yscale('log')

for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(18)

rect = Rectangle((2.5,1e-4),17,1e4,linewidth=1,edgecolor='none',facecolor='grey', alpha = 0.7)
ax1.add_patch(rect)
plt.vlines(2.5, 1e-4, 1, colors='r')
#plt.show()
plt.savefig('Sobol2T.pdf', bbox_inches='tight') #DEBUG


#Keeping only 4 vars
fig, ax1 = plt.subplots(1,1)
ax1.scatter(range(majI), sVal[:majI], color='k', marker='o', facecolors='none')
#ax1.errorbar(range(34), firstOrder, yerr=firstError, fmt='x', ecolor='k', ls='none') #, marker='+', label=r'first order indices')
#ax1.errorbar(range(34), totalOrder, yerr=totalError, ecolor='r') #, marker='o', facecolors='none', label=r'total order indices')

ax1.set_ylim([1e-4, 1e0])
ax1.set_xticks(range(majI))
ax1.set_xticklabels(sKey[:majI], rotation='vertical', fontsize=18)

ax1.set_ylabel(r'Sobol index')
plt.yscale('log')

for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
  item.set_fontsize(18)

rect = Rectangle((3.5,1e-4),16,1e4,linewidth=1,edgecolor='none',facecolor='grey', alpha = 0.7)
ax1.add_patch(rect)
plt.vlines(3.5, 1e-4, 1, colors='r')
#plt.show()
plt.savefig('Sobol3.pdf', bbox_inches='tight') #DEBUG













#------------------ END TEST

prt("THE END", 'green', True)
quit()
