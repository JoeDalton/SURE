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
parser.add_argument('--noWrite',    action='store_false',       help = "If present, write to file")
args        = parser.parse_args()
shouldWrite = args.noWrite

# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'SobolIndicesFinder.input'
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
  DumpObject(toDump, outFile, outPath)
else:
  print(firstOrder)



prt("THE END", 'green', True)
quit()
