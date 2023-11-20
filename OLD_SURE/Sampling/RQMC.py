#from Utils.I_O import prt
from Sampling  import Sampling
from QMC       import QMC
import openturns as ot
import time
#import pickle


#TODO: Remove this class as it may be handled directly in the CV Run Controller ?


class RQMC(Sampling):
  
  ######################
  #     Properties     #
  ######################
  sequence           = ''
  firstSample        = 0
  nSample            = 0
  nQMCDrawer         = 0
  xQMCDrawer         = []

  ###################
  #     Methods     #
  ###################
  def __init__(self, parent, distribution, nSample, nQMCDrawer=10, firstSample=0, sequence='Sobol', verbosity=-1, fixSeed=False):
    self.parent       = parent
    self.verbosity    = verbosity
    self.nSample      = nSample
    self.firstSample  = firstSample
    self.sequence     = sequence
    
    if fixSeed:
      ot.RandomGenerator.SetSeed(0)
    else:
      ot.RandomGenerator.SetSeed(int(1000*time.time()))

    for i in range(nQMCDrawer):
      d = QMC(self, distribution, nSample, firstSample, sequence, verbosity=self.verbosity-1)
      d.setRandomize(True)
      self.xQMCDrawer.append(d)
  #def Load(self):
  #  pass

  #def Copy(self, other):
  #  pass

  def Draw(self):
    #WARNING: HERE, the "samples" output variable is a list of "ot.samples" objects, while other Sampling drawers output one "ot.samples" object
    samples = []
    for d in self.xQMCDrawer:
      samples.append(d.Draw())
    
    return samples
  
