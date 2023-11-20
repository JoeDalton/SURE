#from Utils.I_O import prt
from Sampling  import Sampling
import openturns as ot
#import pickle
import time
import numpy as np

class MC(Sampling):
  
  ######################
  #     Properties     #
  ######################
  nSample            = 0
  firstSample        = 0
  drawer             = None

  ###################
  #     Methods     #
  ###################
  def __init__(self, parent, distribution, nSample, firstSample=0, verbosity=-1, fixSeed=False):
    if fixSeed:
      ot.RandomGenerator.SetSeed(0)
    else:
      ot.RandomGenerator.SetSeed(int(1000*time.time()))
    self.parent       = parent
    self.verbosity    = verbosity
    self.nSample      = nSample
    self.firstSample  = firstSample
    self.drawer       = ot.MonteCarloExperiment(distribution, int(nSample+firstSample))

  #def Load(self):
  #  pass

  #def Copy(self, other):
  #  pass

  def Draw(self):
    tempSamples = self.drawer.generate()
    
    if self.firstSample==0:
      self.xPoint = np.array(tempSamples)
      return tempSamples
    else:
      #TODO: See how to handle this case properly
      #if dimension==1:
      #  prt('Sampling/Sampling.py/Draw: 1D experiments have not been tested. Please try again later', 'red', self.verbosity)
      #  return None
      #else:
      samples = tempSamples[self.firstSample:,:]
      self.xPoint = np.array(samples)
      return samples
  
