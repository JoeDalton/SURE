import h5py
import numpy as np
import sys
from Utils.I_O import prt
from Utils.I_O import DictFromH5
from Interface.Modifications import *
from VarTransformation import *

class Variable:
    
  ######################
  #     Properties     #
  ######################
  ID                  = ''
  savedItems          = {}
  parent              = None
  verbosity           = -1
  ###################
  #     Methods     #
  ###################
  def __init__(self):
    pass

  def Load(self):
  	pass
  

  def Copy(self, other):
  	pass
  

  def Dump(self, directory='./'):
    DumpObject(self.savedItems, directory)
    prt(self.redType + "_" + self.ID + " : Dumped", "green", self.verbosity)

    

class AbstractVariable(Variable):

  ######################
  #     Properties     #
  ######################
  ###################
  #     Methods     #
  ###################
  def __init__(self):
  	pass
  

  def Load(self):
  	pass


  def Save(self):
  	pass



class PhysicalVariable(Variable):

  ######################
  #     Properties     #
  ######################
  targetFile            = ''
  targetDictPath        = ''
  xTransformationParam  = []
  transformationType    = ''
  Tranformation         = None
  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    self.savedItems = {}
    if 'empty' in kwargs and 'verbosity' in kwargs:
      self.verbosity              = kwargs['verbosity']
    else:
      willExit = False
      try:
        self.ID                   = kwargs['ID']                    # String - Variable name
        self.targetFile           = kwargs['targetFile']            # String - Relative path of file
        self.targetDictPath       = kwargs['targetDictPath']        # String - Path in the target file, with / as separator
        self.transformationType   = kwargs['transformationType']    # String - type of variable physical transformation
        self.xTransformationParam = kwargs['xTransformationParam']  # List of float - parameters of physical transformation
        self.verbosity            = kwargs['verbosity']             # integer - verbosity level
      except:
        prt('', 'red', True)
        prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
        prt('Error: The arguments must contain either ("empty") or ("ID", "parent", "tranformationType", "xTranformationParam", "targetFile", "targetDictPath", "verbosity").', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit()
        
      if self.transformationType == 'NormalToLogNormal':
        # xTransformationParam must be: [nominalValue, uncertaintyFactor]
        if len(self.xTransformationParam) == 2:
          self.Transformation = NormalToLogNormal
        else:
          prt('', 'red', True)
          prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
          prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape [nominalValue, uncertaintyFactor]', 'red', True)
          prt('Exiting', 'red', True)
          sys.exit()
      elif self.transformationType == 'StretchedUniform':
        # xTransformationParam must be: [minValue, maxValue]
        if len(self.xTransformationParam) == 2:
          self.Transformation = StretchedUniform
        else:
          prt('', 'red', True)
          prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
          prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape [minValue, maxValue]', 'red', True)
          prt('Exiting', 'red', True)
          sys.exit()
      elif self.transformationType == 'UniformToLogUniform':
        # xTransformationParam must be: [minValue, maxValue]
        if len(self.xTransformationParam) != 2:
          prt('', 'red', True)
          prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
          prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape [minValue, maxValue]', 'red', True)
          prt('Exiting', 'red', True)
          sys.exit()
        elif self.xTransformationParam[0] == 0.0:
          prt('', 'red', True)
          prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
          prt('Error: For ' + self.transformationType + ', minValue can not be 0 !', 'red', True)
          prt('Exiting', 'red', True)
          sys.exit()
        else:
          self.Transformation = UniformToLogUniform
      elif self.transformationType == 'None':
        # xTransformationParam must be: []
        if len(self.xTransformationParam) == 0:
          self.Transformation = NoTransform
        else:
          prt('', 'red', True)
          prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
          prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape []', 'red', True)
          prt('Exiting', 'red', True)
          sys.exit()
      elif self.transformationType == 'Dirac':
        # xTransformationParam must be: [Value]
        if len(self.xTransformationParam) == 1:
          self.Transformation = Dirac
        else:
          prt('', 'red', True)
          prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
          prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape [value]', 'red', True)
          prt('Exiting', 'red', True)
          sys.exit()
      else:
        prt('', 'red', True)
        prt('Experiment/Variable.py/PhysicalVariable/__init__', 'red', True)
        prt('Error: Invalid argument "transformationType": ' + self.transformationType + '. Valid arguments are "Dirac", "NormalToLogNormal", "StretchedUniform", and "None"', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit()

  
  def ModifyTarget(self, value, workdir):
    targetValue = self.Transformation(value, self.xTransformationParam)
    if type(self.targetFile) is list:
      if type(self.targetDictPath) is list:
        if len(self.targetFile) == len(self.targetDictPath):
          for i in range(len(self.targetDictPath)):
            ModifyTargetFile(workdir + '/' + self.targetFile[i], self.targetDictPath[i], targetValue)
        else:
          prt('Error, inconsistency between targetFile and targetDictPath 1', 'red', True)
          print(self.targetFile)
          print(self.targetDictPath)
          sys.exit()
      else:
        prt('Error, inconsistency between targetFile and targetDictPath 2', 'red', True)
        sys.exit()
    else:
      if type(self.targetDictPath) is list:
        prt('Error, inconsistency between targetFile and targetDictPath 3', 'red', True)
        sys.exit()
      else: 
        ModifyTargetFile(workdir + '/' + self.targetFile, self.targetDictPath, targetValue)


  def Save(self):
    self.savedItems['ID']                   = self.ID
    self.savedItems['ObjectType']           = 'PhysicalVariable'
    self.savedItems['xTransformationParam'] = self.xTransformationParam
    self.savedItems['targetFile']           = self.targetFile
    self.savedItems['targetDictPath']       = self.targetDictPath
    self.savedItems['transformationType']   = self.transformationType


  def Load(self, myInput, prefix=''):
    if isinstance(myInput, dict):
      self._LoadFromDictionary(myInput)
    elif isinstance(myInput, str):
      self._LoadFromFile(myInput, prefix=prefix)
    else:
      prt('', 'red', True)
      prt('Experiment/Variable.py.py/PhysicalVariable/Load', 'red', True)
      prt('The argument "myInput" must be of type "dict" or "str" (name of file)', 'red', True)
      prt('Exiting...', 'red', True)
      sys.exit()


  def _LoadFromFile(self, myFile, prefix=''):
    prt('', 'green', self.verbosity)
    prt('Loading the variable from the input file ' + myFile +'.', 'green', self.verbosity)
    myDict = DictFromH5(myFile, prefix=prefix)
    self._LoadFromDictionary(myDict)


  def _LoadFromDictionary(self, myDict):
    if myDict['ObjectType'] != 'PhysicalVariable':
      prt('', 'red', True)
      prt('Experiment/Variable.py/PhysicalVariable/_LoadFromDictionnary', 'red', True)
      prt('The provided dictionnary defines an object of type ' + str(myDict['ObjectType']) + ' while it should define an object of type PhysicalVariable', 'red', True)
      prt('Exiting...', 'red', True)
      sys.exit()
    self.ID                   = myDict['ID']
    self.xTransformationParam = myDict['xTransformationParam']
    self.targetFile           = myDict['targetFile']
    self.targetDictPath       = myDict['targetDictPath']
    self.transformationType   = myDict['transformationType']

    if self.transformationType == 'NormalToLogNormal':
      # xTransformationParam must be: [nominalValue, uncertaintyFactor]
      if len(self.xTransformationParam) == 2:
        self.Transformation = NormalToLogNormal
      else:
        prt('', 'red', True)
        prt('Experiment/Variable.py/PhysicalVariable/_LoadFromDictionary', 'red', True)
        prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape [nominalValue, uncertaintyFactor]', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit()
    elif self.transformationType == 'StretchedUniform':
      # xTransformationParam must be: [minValue, maxValue]
      if len(self.xTransformationParam) == 2:
        self.Transformation = StretchedUniform
      else:
        prt('', 'red', True)
        prt('Experiment/Variable.py/PhysicalVariable/_LoadFromDictionary', 'red', True)
        prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape [minValue, maxValue]', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit()
    elif self.transformationType == 'None':
      # xTransformationParam must be: []
      if len(self.xTransformationParam) == 0:
        self.Transformation = NoTransform
      else:
        prt('', 'red', True)
        prt('Experiment/Variable.py/PhysicalVariable/_LoadFromDictionary', 'red', True)
        prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape []', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit()
    elif self.transformationType == 'Dirac':
      # xTransformationParam must be: [Value]
      if len(self.xTransformationParam) == 1:
        self.Transformation = Dirac
      else:
        prt('', 'red', True)
        prt('Experiment/Variable.py/PhysicalVariable/_LoadFromDictionary', 'red', True)
        prt('Error: For ' + self.transformationType + ', xTransformationParam must be of shape [Value]', 'red', True)
        prt('Exiting', 'red', True)
        sys.exit()
    else:
      prt('', 'red', True)
      prt('Experiment/Variable.py/PhysicalVariable/_LoadFromDictionary', 'red', True)
      prt('Error: Invalid argument "transformationType": ' + self.transformationType + '. Valid arguments are "Dirac", "NormalToLogNormal", "StretchedUniform", and "None"', 'red', True)
      prt('Exiting', 'red', True)
      sys.exit()

    prt('Variable ' + self.ID + ' successfully loaded.', 'green', self.verbosity)
