#!/usr/bin/env python2.7
import os
import sys
import h5py
import Link
import argparse
import openturns                                  as ot
from Utils.I_O                import prt
from Utils.I_O                import InputHandler as IH
from Experiment.Experiment    import Experiment
from Experiment.RunController import StdRunCtrl
from Interface.Simulation     import Simulation

# ----- Parsing arguments -----
parser = argparse.ArgumentParser(description="RunExperiment")
parser.add_argument('-i',	  required=False,	type=str,  help = "Relative path of input file")
args = parser.parse_args()


# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'RunExperiment.input'
  try:
    INPUT = IH(inputFile, modify=False)
  except:
    prt('Input File does not exist', 'None', True)
    prt('  --> Copy an example in ' + inputFile, 'None', True)
    os.system("cp " + Link.sureHome + '/InputFiles/' + inputFile + ' .')
    sys.exit()

# Simulation-related reading
INPUT.GetBlock('Simulation')
mySimulation  = Simulation(block=INPUT)

# Experiment-related reading
INPUT.GetBlock('Experiment')
expPath       = INPUT.GetCharacterFromKey('ExpPath')
myExperiment  = Experiment(block=INPUT, simulation=mySimulation)

# Controller-related reading
INPUT.GetBlock('Controller')
myController = StdRunCtrl(block=INPUT, experiment=myExperiment)

# ----- Running controller -----
myController.PrepareSamples()
myController.Run()

# ----- Saving the experiment object ----
if not myController.isContinued:
  myExperiment.Save()
  myExperiment.Dump(fileName=myController.resultFile, path=expPath)

prt("THE END", 'green', True)
quit()
