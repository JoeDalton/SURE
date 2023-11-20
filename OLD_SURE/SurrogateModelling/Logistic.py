import numpy as np
import pickle
import sys
from sklearn.linear_model import LogisticRegression
from Utils.I_O import prt, ErrorMsg
from Utils.I_O                        import DictFromH5
from SurrogateModelling.SurrogateModel import SurrogateModel as SM

fileName = 'SurrogateModelling/Logistic.py'

class Logistic(SM):
  
  ######################
  #     Properties     #
  ######################
  responseSurface=None
  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    self.savedItems   = {}
    willExit = True
    self.objectType   = 'Surrogate_Logistic'
    if "block" in kwargs and "experiment" in kwargs: # kwargs["block"] must be a block of "InputHandler"
      self.parent       = kwargs["experiment"] # Must be a valid experiment object with valid variables and distribution
      self._InitWithBlock(kwargs["block"])
      willExit          = False
    elif "empty" in kwargs and "verbosity" in kwargs:
      self.verbosity    = kwargs["verbosity"]
      willExit          = False
    else:
      try:
        self.ID           = kwargs["ID"]
        self.verbosity    = kwargs["verbosity"]
        if "experiment" in kwargs:
          self.parent       = kwargs["experiment"]
        willExit          = False
      except:
        willExit = True
    if willExit:
      ErrorMsg('Error: The arguments must contain ("block" and "experiment") or ("ID", "experiment" and "verbosity") or ("ID" and "verbosity") or ("empty=True" and "verbosity").', fileName)


  def _InitWithBlock(self, BLOCK):
    prt('', 'green', self.verbosity)
    prt('Initializing the Logistic regression from the input file.', 'green', self.verbosity)
    self.ID         = BLOCK.GetCharacterFromKey('ID')
    self.verbosity  = BLOCK.GetCharacterFromKey('outputLevel')


  def Load(self, myInput, prefix=''):
    if isinstance(myInput, dict):
      self._LoadFromDictionary(myInput)
    elif isinstance(myInput, str):
      self._LoadFromFile(myInput, prefix=prefix)
    else:
      ErrorMsg('The argument "myInput" must be of type "dict" or "str" (name of file)', fileName)


  def _LoadFromFile(self, myFile, prefix=''):
    prt('', 'green', self.verbosity)
    prt('Loading the Logistic regression from the input file ' + myFile +'.', 'green', self.verbosity)
    myDict = DictFromH5(myFile, prefix=prefix)
    self._LoadFromDictionary(myDict)


  def _LoadFromDictionary(self, myDict):
    if myDict['ObjectType'] != self.objectType:
      ErrorMsg('The provided dictionary defines an object of type ' + str(myDict['ObjectType']) + ' while it should define an object of type ' + self.objectType, fileName)
    self.ID = myDict['ID']
    self.responseSurface  = pickle.loads(myDict['ResponseSurface'])
    prt('Logistic regression ' + self.ID + ' successfully loaded.', 'green', self.verbosity)


  def Save(self):
    self.savedItems['ID']              = self.ID
    self.savedItems['ObjectType']      = self.objectType
    self.savedItems['ResponseSurface'] = pickle.dumps(self.responseSurface)


  def Build(self, xSample, xEval, debug=False):
    if debug:
      self.responseSurface = LogisticRegression(random_state=0)
    else:
      self.responseSurface = LogisticRegression()
    self.responseSurface.fit(xSample, xEval)


  def Evaluate(self, inputs, proba=False):
    if proba:
      return self.responseSurface.predict_proba(inputs)
    else:
      return self.responseSurface.predict(inputs)