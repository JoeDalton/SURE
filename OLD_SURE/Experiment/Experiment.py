import h5py
import numpy        as np
import sys
import pickle
import openturns    as ot
from Utils.I_O            import prt
from Utils.I_O            import DumpObject
from Utils.I_O            import DictFromH5
from Variable             import *
from Interface.Simulation import Simulation

class Experiment:
    
  ######################
  #     Properties     #
  ######################
  ID                 = ''
  savedItems         = {}
  xVariable          = []
  xVarID             = []
  verbosity          = -1
  simulation         = None
  distribution       = None
  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    willExit = True
    self.savedItems   = {}
    if "block" in kwargs and "simulation" in kwargs: # kwargs["block"] must be a block of "InputHandler"
      if "simulation" in kwargs:
        self.simulation   = kwargs["simulation"]
      else:
        self.simulation   = None
      self._InitWithBlock(kwargs["block"]) 
      willExit            = False
    elif "empty" in kwargs and "verbosity" in kwargs:
      self.verbosity      = kwargs['verbosity']
      willExit            = False
    else:
      try:
        self.ID           = kwargs["ID"]
        self.verbosity    = kwargs["verbosity"]
        self.simulation   = kwargs["simulation"]
        self.xVariable    = []          # List of variable objects
        self.distribution = None
        willExit          = False
      except:
        willExit          = True

    if willExit:
      prt('', 'red', True)
      prt('Experiment/Experiment.py/Experiment/__init__', 'red', True)
      prt('Error: The arguments must contain either ("block") or ("block" and "simulation") or ("empty" and "verbosity") or ("ID", "simulation" and "verbosity").', 'red', True)
      prt('Exiting', 'red', True)
      sys.exit() 


  def Copy(self, other):
    pass
 
 
  def Dump(self, fileName="dummy.h5", path=''):
    DumpObject(self.savedItems, fileName, path)
    prt("Experiment_" + self.ID + " : Dumped", "green", True) #, self.verbosity)
    

  def SetVariables(self, xVar, distribution):
    self.xVariable    = xVar
    self.distribution = distribution


  def _InitWithBlock(self, BLOCK):
    keys = BLOCK.GetKeysForBlock()
    if 'LoadFile' in keys:
      myFile = BLOCK.GetCharacterFromKey('LoadFile')
      if 'LoadPath' in keys:
        myPath = BLOCK.GetCharacterFromKey('LoadPath')
      else:
        myPath = ''
      self.Load(myFile, prefix=myPath)
    else:
      self.ID             = BLOCK.GetCharacterFromKey('ID')
      self.xVarID         = []
      xVar                = []
      xTargetFile         = []
      xTargetDictPath     = []
      xDistribution       = []
      xTransformationType = []
      xxTransformParam    = []

      # Variable-related reading
      tempList      = BLOCK.GetSubBlockNames()
      metaVarList   = []
      for item in tempList:
        if 'Variable' in item: # --------------------------SUBBLOCKs containing variables must have 'Variable' in their name
          metaVarList.append(item)

      for metaVar in metaVarList:
        BLOCK.GetSubBlock(metaVar)

        varKeys = BLOCK.GetKeysForBlock()
        mode              = BLOCK.GetCharacterFromKey('mode')
        if mode == 'Multiple':
          if 'idxList' in varKeys:
            xVarNum       = BLOCK.GetIntegerFromKey('idxList')
            if 'number' in varKeys or 'firstIdx' in varKeys:
              prt('', 'yellow', self.verbosity)
              prt('Experiment/Experiment.py/Experiment/_InitWithBlock:', 'yellow', self.verbosity)
              prt('Warning: Keywords "number" and "firstIdx" are ignored in the presence of "idxList"', 'yellow', self.verbosity)
          elif 'number' in varKeys:
            number        = BLOCK.GetIntegerFromKey('number')
            if 'firstIdx' in varKeys:
              firstIdx    = BLOCK.GetIntegerFromKey('firstIdx')
            else:
              firstIdx    = 0
            xVarNum       = range(firstIdx, number + firstIdx)
          else:
            prt('', 'red', True)
            prt('Experiment/Experiment.py/Experiment/_InitWithBlock', 'red', True)
            prt('Error: Blocks of multiple variables must contain either "idxList" or "number"', 'red', True)
            prt('Exiting', 'red', True)
            sys.exit()
        elif mode == 'Single':
          xVarNum         = [0]
          firstIdx    = 0
        else:
          prt('', 'red', True)
          prt('Experiment/Experiment.py/Experiment/_InitWithBlock', 'red', True)
          prt('Error: Variable mode is either "Multiple" or "Single".', 'red', True)
          prt('Exiting', 'red', True)
          sys.exit()

        for varIdx in xVarNum:
          self.xVarID.append(BLOCK.GetCharacterFromKey('ID').replace('%',str(varIdx)))

          # Special treatment for the case where there are several files in which the variable must be modified
          # TODO: Check here the consistency between the length (where appropriate) of the lists of target files and target paths
          tempTargetFile      = BLOCK.GetCharacterFromKey('targetFile')
          tempTargetDictPath  = BLOCK.GetCharacterFromKey('targetDictPath')
          newTargetFile       = []
          newTargetDictPath   = []
          if type(tempTargetFile) is list:
            for item in tempTargetFile:
              newTargetFile.append(item.replace('%',str(varIdx)))
            xTargetFile.append(newTargetFile)
          else:
            xTargetFile.append(tempTargetFile.replace('%',str(varIdx)))
          if type(tempTargetDictPath) is list:
            for item in tempTargetDictPath:
              newTargetDictPath.append(item.replace('%',str(varIdx)))
            xTargetDictPath.append(newTargetDictPath)
          else:
            xTargetDictPath.append(tempTargetDictPath.replace('%',str(varIdx)))
          

          xDistribution.append(BLOCK.GetCharacterFromKey('distribution').replace('%',str(varIdx)))
          xTransformationType.append(BLOCK.GetCharacterFromKey('transformationType').replace('%',str(varIdx)))

          transformFrom = BLOCK.GetCharacterFromKey('transformationFrom')
          if transformFrom == 'Values':
            if xTransformationType[varIdx-firstIdx] == 'Dirac':
              tempVal = BLOCK.GetRealFromKey('transformationValues')
              if mode == 'Single':
                xxTransformParam.append([tempVal])
              else:
                xxTransformParam.append([tempVal[varIdx-firstIdx]])
            elif xTransformationType[varIdx-firstIdx] == 'None':
              xxTransformParam.append([])
            else:
              xxTransformParam.append(BLOCK.GetRealFromKey('transformationValues'))
              # TODO : This is not well defined for several variables with given transformationValues !
          elif transformFrom == 'Files':
            temp = []
            xFile = BLOCK.GetCharacterFromKey('transformationFiles')
            xPath = BLOCK.GetCharacterFromKey('transformationPaths')
            for i in range(len(xFile)):
              xFile[i] = xFile[i].replace('%',str(varIdx))
              xPath[i] = xPath[i].replace('%',str(varIdx))
            #TODO: Clean this to handle any type of file + check that xFile and xPath have the same length
            for fileIdx in range(len(xFile)):
              myFile = xFile[fileIdx]
              myPath = xPath[fileIdx]
              f = h5py.File(myFile, 'r')
              try:
                temp.append(f[myPath][0])
              except ValueError:
                temp.append(f[myPath][()])
              f.close()
            xxTransformParam.append(temp)
          else:
            prt('Invalid "transformFrom" argument', 'red', True) #TODO: better error handling

      # ----- Initializing joint distribution -----
      distributionList = []
      for distrib in xDistribution:
        if distrib == 'Normal':
          distributionList.append(ot.Normal(0.0,1.0))
        elif distrib == 'Uniform':
          distributionList.append(ot.Uniform(0.0,1.0))
        else:
          prt('Unknown distribution type (' + distrib + '). Supported types are "Normal" and "Uniform".', 'red', True)
          prt('Exiting...', 'red', True)
          sys.exit()
      distribution = ot.ComposedDistribution(distributionList)

      # ----- Initializing variables -----
      for varIdx in range(len(self.xVarID)):
        var = PhysicalVariable(ID=self.xVarID[varIdx], parent=self, \
                                transformationType=xTransformationType[varIdx], \
                                xTransformationParam=xxTransformParam[varIdx], \
                                targetFile=xTargetFile[varIdx], \
                                targetDictPath=xTargetDictPath[varIdx], \
                                verbosity=self.verbosity-1)
        xVar.append(var)

      # ----- Setting the variables of the experiment -----
      self.SetVariables(xVar, distribution)


  def Save(self):
    # Make sure all the variables are saved and build a dictionary of saved variables
    for var in self.xVariable:
      var.Save()
    varDict = {}
    for vIdx in range(len(self.xVariable)):
      varDict[self.xVarID[vIdx]] = self.xVariable[vIdx].savedItems
    # Make sure the underlying simulation object is saved
    self.simulation.Save()
    # Save self
    self.savedItems['ID']                   = self.ID
    self.savedItems['ObjectType']           = 'Experiment'
    self.savedItems['Variables']            = varDict
    self.savedItems['Variable_IDs']         = self.xVarID # Necessary to keep the right order of variables
    self.savedItems['Distribution']         = pickle.dumps(self.distribution)
    self.savedItems['Simulation']           = self.simulation.savedItems


  def Load(self, myInput, prefix=''):
    if isinstance(myInput, dict):
      self._LoadFromDictionary(myInput)
    elif isinstance(myInput, str):
      self._LoadFromFile(myInput, prefix=prefix)
    else:
      prt('', 'red', True)
      prt('Experiment/Experiment.py/Experiment/Load', 'red', True)
      prt('The argument "myInput" must be of type "dict" or "str" (name of file)', 'red', True)
      prt('Exiting...', 'red', True)
      sys.exit()


  def _LoadFromFile(self, myFile, prefix=''):
    prt('', 'green', self.verbosity)
    prt('Loading the experiment from the input file ' + myFile +'.', 'green', self.verbosity)
    myDict = DictFromH5(myFile, prefix=prefix)
    self._LoadFromDictionary(myDict)


  def _LoadFromDictionary(self, myDict):
    if myDict['ObjectType'] != 'Experiment':
      prt('', 'red', True)
      prt('Experiment/Experiment.py/Experiment/_LoadFromDictionnary', 'red', True)
      prt('The provided dictionnary defines an object of type ' + str(myDict['ObjectType']) + ' while it should define an object of type Experiment', 'red', True)
      prt('Exiting...', 'red', True)
      sys.exit()
    # Load self
    self.ID           = myDict['ID']
    self.distribution = pickle.loads(myDict['Distribution'])
    self.xVarID       = myDict['Variable_IDs']
    
    # Load variables
    self.xVariable = []
    for varID in self.xVarID:
      if myDict['Variables'][varID]['ObjectType'] == 'PhysicalVariable':
        var = PhysicalVariable(empty=True, verbosity=self.verbosity-1)
      else:
        prt('', 'red', True)
        prt('Experiment/Experiment.py/Experiment/_LoadFromDictionnary', 'red', True)
        prt('Variable ' + varID + ': The object type ' + myDict['Variables'][varID]['ObjectType'] + ' is not yet supported.', 'red', True)
        prt('Exiting...', 'red', True)
        sys.exit()
      var.Load(myDict['Variables'][varID])
      self.xVariable.append(var)

        
    # Load simulation
    self.simulation = Simulation(empty=True, verbosity=self.verbosity-1)
    self.simulation.Load(myDict['Simulation'])
 
    prt('Experiment ' + self.ID + ' successfully loaded.', 'green', self.verbosity)
