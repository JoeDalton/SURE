#from Utils.I_O import prt
from Utils.I_O import DumpObject

class Sampling:
  
  ######################
  #     Properties     #
  ######################
  parent             = None
  SavedItems         = {}
  verbosity          = -1
  xPoint             = None # np.array
  #xWeight            = None # np.array Used only for deterministic samplings TODO:clean ?
  ###################
  #     Methods     #
  ###################
  def __init__(self, parent=None,verbosity=-1):
    self.parent     = parent
    self.verbosity  = verbosity

  def Dump(self, directory="./"):
    DumpObject(self.SavedItems, directory)
    prt(self.objectType + "_" + self.ID + " : Dumped", "green", self.verbosity)

