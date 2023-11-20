from Utils.I_O import prt, ErrorMsg
from Utils.I_O import DumpObject

fileName = "SurrogateModelling/SurrogateModel.py"

class SurrogateModel:
  
  ######################
  #     Properties     #
  ######################
  parent             = None # Experiment object
  savedItems         = {}
  objectType         = ""
  ID                 = ""
  verbosity          = False
  ###################
  #     Methods     #
  ###################
  #def __init__(self, parent, ID):
  #  pass
  #def Load(self, dataFile):
  #  pass

  def Copy(self, other):
    pass

  def Evaluate(self, inputs):
    ErrorMsg('This surrogate model has no "Evaluate" function', fileName)
    return None

  def Dump(self, fileName="dummy.h5", path=''):
    DumpObject(self.savedItems, fileName, path)
    prt(self.objectType + "_" + self.ID + " : Dumped", "green", self.verbosity)

  def __call__(self, sample):
    try:
      return self.Evaluate(sample)[0]
    except:
      return self.Evaluate(sample)
      
