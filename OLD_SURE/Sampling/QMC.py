from Utils.I_O import prt, ErrorMsg
from Sampling  import Sampling
import openturns as ot
#import pickle
import numpy as np
import sys

fileName = 'Sampling/QMC.py'

class QMC(Sampling):
  
  ######################
  #     Properties     #
  ######################
  sequence           = ''
  firstSample        = 0
  nSample            = 0
  drawer             = None

  ###################
  #     Methods     #
  ###################
  def __init__(self, parent, distribution, nSample, firstSample=0, sequence='Sobol', verbosity=-1):
    self.parent       = parent
    self.verbosity    = verbosity
    self.nSample      = nSample
    self.firstSample  = firstSample
    self.sequence     = sequence

    if self.sequence == 'LHS':
      self.drawer = ot.LHSExperiment(distribution, int(nSample+firstSample))
    else:
      if self.sequence == 'Sobol':
        seq = ot.SobolSequence()
      elif self.sequence == 'Halton':
        seq = ot.HaltonSequence()
      elif self.sequence == 'Reverse Halton':
        seq = ot.ReverseHaltonSequence()
      elif self.sequence == 'Faure':
        seq = ot.FaureSequence()
      elif self.sequence == 'Haselgrove':
        seq = ot.HaselgroveSequence()
      else:
        ErrorMsg('Error: Please provide a valid "sequence" argument. Valid arguments are "Sobol", "Halton", "Reverse Halton", "Faure", "Haselgrove" and "LHS"', fileName)
        return None
      self.drawer = ot.LowDiscrepancyExperiment(seq, distribution, int(self.nSample+self.firstSample), False)

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
  
  def setRandomize(self, isRandom):
    self.drawer.setRandomize(isRandom)
