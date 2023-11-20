#!/usr/bin/env python2.7

import os
import sys
import h5py
import Link
import argparse
import numpy                                                      as np
import openturns                                                  as ot
import matplotlib.pyplot                                          as plt
import Utils.PrettyPlots                                          as pp
from Utils.I_O                        import prt, DumpObject, GetObjectType
from Utils.I_O                        import InputHandler         as IH
from Sampling.QMC                     import QMC
from Experiment.Experiment            import Experiment
from SurrogateModelling.PCE           import *
from SurrogateModelling.OneDSpline    import *
from SurrogateModelling.Kriging       import *
from DimensionReduction.Preprocessing import Preprocessor


# ----- Parsing arguments -----
parser      = argparse.ArgumentParser(description="SurrogateResampler")
parser.add_argument('-i',	          required=False,	type=str,   help = "Relative path of input file")
args        = parser.parse_args()

# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'SurrogateResampler.input'
  try:
    INPUT = IH(inputFile, modify=False)
  except:
    prt('Input File does not exist', 'None', True)
    prt('  --> Copy an example in ' + inputFile, 'None', True)
    os.system("cp " + Link.sureHome + '/InputFiles/' + inputFile + ' .')
    sys.exit()



#TODO: Add the preprocessing and postprocessing of samples !!!

##############################
#     INPUT FILE READING     #
##############################
# ----- Surrogate-related inputs -----
INPUT.GetBlock('Surrogate')
keys            = INPUT.GetKeysForBlock()
surFile         = INPUT.GetCharacterFromKey('SurrogateFile')
surPath         = INPUT.GetCharacterFromKey('SurrogatePath')
ExpeFile        = INPUT.GetCharacterFromKey('ExperimentFile')
ExpePath        = INPUT.GetCharacterFromKey('ExperimentPath')
surType         = GetObjectType(surFile, surPath)

#try:
#  INPUT.GetBlock('Preprocessing')
#  isPreProc = True
#except:
#  isPreproc = False
isPreProc = False #DEBUG

##############################
#        DATA LOADING        #
##############################
# ----- Loading experiment -----
myExp = Experiment(empty=True, verbosity=1)
myExp.Load(ExpeFile, prefix=ExpePath)

# ----- Loading surrogate -----
if surType == 'Surrogate_PCE':
  mySur = PCE(empty=True, verbosity=1)
elif surType == 'Surrogate_Kriging':
  mySur = Kriging(empty=True, verbosity=1)
elif surType == 'Surrogate_1DSpline':
  mySur = OneDSpline(empty=True, verbosity=1) 
else:
  print('Surrogate Type "' + surType + '" not implemented.')
  quit()
mySur.Load(surFile, prefix=surPath)

distribution  = myExp.distribution
nSample       = 100000
xEval         = np.zeros(nSample)
MCDrawer      = QMC(None, distribution, nSample, firstSample=0, verbosity=1, sequence='LHS')
MCDrawer.Draw()
dim           = distribution.getDimension()
#SanityCheck
MCDrawer.xPoint[:,2] = MCDrawer.xPoint[0,2]
xSample       = MCDrawer.xPoint

if isPreProc:
  INPUT.GetBlock('Preprocessing')
  xSample, dummy = Preprocessor(xSample, None, INPUT)

if mySur.objectType == 'Surrogate_1DSpline':
  xEval = mySur.Evaluate(xSample)
else:
  for i in range(nSample):
    xEval[i] = mySur.Evaluate(xSample[i])[0]

if isPreProc:
  INPUT.GetBlock('Preprocessing')
  dummy, xEval = Preprocessor(None, xEval, INPUT)


#resampPath = surPath + "/Resampling"
#SanityCheck
resampPath = surPath + "/Sanity_Check"
dict={}
dict['Data'] = xEval
DumpObject(dict, fileName=surFile, prefix=resampPath)



prt("THE END", 'green', True)
quit()
