from Utils.I_O import prt
from Utils.I_O import DumpObject

class DimensionReduction:
  
  ######################
  #     Properties     #
  ######################
  parent             = None
  savedItems         = {}
  verbosity          = -1
  totalDim           = 0
  reducedDim         = 0
  objectType         = ''
  ID                 = ''

  ###################
  #     Methods     #
  ###################
  def __init__(self, ID, parent=None, verbosity=-1):
    self.parent     = parent
    self.verbosity  = verbosity
    self.savedItems = {}
    self.totalDim   = 0
    self.reducedDim = 0
    self.objectType = ''
    self.ID         = ID

  def Dump(self, fileName="dummy.h5", path=''):
    DumpObject(self.savedItems, fileName, path)
    prt(self.objectType + "_" + self.ID + " : Dumped", "green", self.verbosity)
