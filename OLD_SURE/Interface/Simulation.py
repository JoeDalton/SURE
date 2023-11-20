import h5py
import numpy as np
import sys
from Utils.I_O import prt, ErrorMsg
from Utils.I_O import DictFromH5, DumpObject
from Interface.Results import *

fileName = 'Interface/Simulation.py'

class Simulation:
    
  ######################
  #     Properties     #
  ######################
  ID                 = ''
  objectType         = ''
  refFolderPath      = ''
  savedItems         = {}
  verbosity          = -1
  RetrieveQoI        = None
  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    self.savedItems = {}
    willExit        = True
    if "block" in kwargs:
      self._InitWithBlock(kwargs["block"]) 
      self.SetRetrieveQoIFunc()    
      willExit        = False
    elif "empty" in kwargs and "verbosity" in kwargs:
      willExit        = False
      self.verbosity  = kwargs['verbosity']
    else:
      try:
        self.ID             = kwargs["ID"]
        self.verbosity      = kwargs["verbosity"]
        self.objectType     = 'Simulation_' + kwargs["simuType"]
        self.refFolderPath  = kwargs["refFolderPath"]
        self.SetRetrieveQoIFunc()    
        willExit        = False
      except:
        willExit = True 
      
      if willExit:
        ErrorMsg('The arguments must contain either ("block") or ("empty" and "verbosity") or ("ID", "simuType", "refFolderPath" and "verbosity").', fileName)



  def SetRetrieveQoIFunc(self):
    if self.objectType == "Simulation_AGATH_0D_IDT":
      self.RetrieveQoI = AGATH_0D_IDT
    elif self.objectType == "Simulation_AGATH_1D_IDT":
      self.RetrieveQoI = AGATH_1D_IDT
    elif self.objectType == "Simulation_To_Keep":
      self.RetrieveQoI = Simulation_To_Keep
    else:
      ErrorMsg('Error: The simulation type ' + self.objectType + ', is not valid or not yet supported. Valid types are "Simulation_AGATH_0D_IDT", "Simulation_AGATH_1D_ID" and "Simulation_To_Keep".', fileName)

  
  def Copy(self, other):
  	pass

  
  def Save(self):
    self.savedItems['ID']         = self.ID
    self.savedItems['ObjectType'] = self.objectType
    self.savedItems['RefFolder']  = self.refFolderPath
 
 
  def Dump(self, directory='./'):
    DumpObject(self.savedItems, directory)
    prt(self.redType + "_" + self.ID + " : Dumped", "green", self.verbosity)
    

  def _InitWithBlock(self, BLOCK):
    self.ID             = BLOCK.GetCharacterFromKey('ID')
    self.objectType     = 'Simulation_' + BLOCK.GetCharacterFromKey('simuType')
    self.refFolderPath  = BLOCK.GetCharacterFromKey('refFolderPath')
    self.verbosity      = -1


  def Load(self, myInput, prefix=''):
    if isinstance(myInput, dict):
      self._LoadFromDictionary(myInput)
    elif isinstance(myInput, str):
      self._LoadFromFile(myInput, prefix=prefix)
    else:
      ErrorMsg('The argument "myInput" must be of type "dict" or "str" (name of file)', fileName)


  def _LoadFromFile(self, myFile, prefix=''):
    prt('', 'green', self.verbosity)
    prt('Loading the simulation from the input file ' + myFile +'.', 'green', self.verbosity)
    myDict = DictFromH5(myFile, prefix=prefix)
    self._LoadFromDictionary(myDict)


  def _LoadFromDictionary(self, myDict):
    #TODO: Something to detect if myDict['ObjectType'] if of form 'Simulation_*'
    #if myDict['ObjectType'] != 'Simulation':
    #  prt('', 'red', True)
    #  prt('Experiment/Experiment.py/Experiment/_LoadFromDictionnary', 'red', True)
    #  prt('The provided dictionnary defines an object of type ' + str(myDict['ObjectType']) + ' while it should define an object of type Experiment', 'red', True)
    #  prt('Exiting...', 'red', True)
    #  sys.exit()
    # Load self
    self.ID             = myDict['ID']
    self.objectType     = myDict['ObjectType']
    self.refFolderPath  = myDict['RefFolder']
    self.SetRetrieveQoIFunc()
    
    prt('Simulation ' + self.ID + ' successfully loaded.', 'green', self.verbosity)
