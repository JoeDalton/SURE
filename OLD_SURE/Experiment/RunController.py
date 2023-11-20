import h5py
import time
import numpy as np
import sys
import os
import subprocess
from Utils.I_O               import prt
from Utils.I_O               import ProgBar
from Interface.Modifications import CopyRunFolder
from Interface.Modifications import EraseRunFolder
from Utils.ResultFileHandler import HandleResultFile
from Sampling.MC             import MC
from Sampling.QMC            import QMC
from Sampling.Deterministic  import CubatureSampling



class RunController:
    
  ######################
  #     Properties     #
  ######################
  parent             = None
  samplerType        = ''
  nThread            = 0
  xPoint             = []
  xEvaluation        = []
  resultFile         = ''
  ID                 = ''
  isContinued        = False
  verbosity          = -1
  ###################
  #     Methods     #
  ###################
  def __init__(self):
  	pass
  
  def Load(self):
  	pass
  
  def Copy(self, other):
  	pass
  
  def Save(self):
  	pass
    


class CVRunCtrl(RunController):
    
  ######################
  #     Properties     #
  ######################
  cvThreshold          = 0
  batchSize            = 0
  nSampler             = 0
  xSampler             = []
  xxPoints             = []
  xxEvaluation         = []
  currentNSample       = 0
  ###################
  #     Methods     #
  ###################
  def __init__(self, ID, parent, samplerType, nSampler, nThread, cvThreshold, batchSize, verbosity):
    # A good compromise for nSampler is 10 for RQMC: Good variance estimation and good cv speed
    self.ID             = ID
    self.parent         = parent
    self.samplerType    = samplerType
    self.nSampler       = nSampler
    self.nThread        = nThread
    self.cvThreshold    = cvThreshold  	
    self.batchSize      = batchSize
    self.currentNSample = 0
    self.xPoint         = []
    self.xEvaluation    = []
    self.verbosity      = verbosity


  def ComputeVariance(self):
    xMeanL = np.zeros(self.nSampler)
    for i in range(self.nSampler):
      xMeanL[i] = np.mean(xxEvaluation[i])
    totMean = np.mean(xMeanL)
    
    var = 0.0
    for i in range(self.nSampler):
      var = (totMean - xMeanL[i]) * (totMean - xMeanL[i])
    var /= (self.nSampler * (self.nSampler - 1.0))
    return var


  def Run(self):
    # Init: Launch the num of runs required by self.batchSize
    # compute var
    # Stock the samples in a list of sample points and a list of evaluation values
    # Loop: 
    # While var > self.cvThreshold:
    #   run a new batch with a StdRunController
    #   compute var
    #   Stock the samples in a list of sample points and a list of evaluation values
    #   
    # TODO: Do it ! Just do it !
    pass









class StdRunCtrl(RunController):
    
  ######################
  #     Properties     #
  ######################
  nSample           = 0
  xLevel            = []
  firstSample       = 0
  sampler           = None
  resultQoIPath     = ''
  resultInputPath   = ''
  resultWeightPath  = ''
  isSparse          = False
  ###################
  #     Methods     #
  ###################
  def __init__(self, **kwargs):
    willExit = True
    self.xPoint           = []
    self.xEvaluation      = []
    if "block" in kwargs and "experiment" in kwargs:
      self.parent       = kwargs["experiment"]
      self._InitWithBlock(kwargs["block"])
      willExit = False
    else:
      #try:
      #  self.ID               = kwargs["ID"]
      #  self.parent           = kwargs["experiment"]
      #  self.nSample          = kwargs["nSample"]
      #  self.resultFile       = kwargs["resultFile"]
      #  self.resultQoIPath    = kwargs["resultQoIPath"]
      #  self.resultInputPath  = kwargs["resultInputPath"]
      #  sequence              = kwargs["sequence"]
      #  self.samplerType      = ''
      #  if sequence == 'None':
      #    self.samplerType  = 'MC'
      #  else:
      #    #TODO: Handle the case of deterministic sampling.
      #    self.samplerType  = 'QMC_' + sequence
      #  if "nThread" in kwargs:
      #    self.nThread      = kwargs["nThread"]
      #  else:
      #    self.nThread      = 1
      #  if "isContinued" in kwargs:
      #    self.isContinued  = kwargs["isContinued"]
      #  else:
      #    self.isContinued  = False
      #  if "verbosity" in kwargs:
      #    self.verbosity    = kwargs["verbosity"]
      #  else:
      #    self.verbosity    = False
      #  willExit = False
      #except:
      #  willExit = True
      #DEBUG TODO: Set this straight
      prt('Direct initialisation of controller is temporarily disabled during development', 'red', True)
      prt('Sorry', 'red', True)
      sys.exit()

    if willExit:
      prt('', 'red', True)
      prt('Experiment/RunController.py/StdRunController/__init__', 'red', True)
      prt('Error: The arguments must contain either ("block" and "experiment") or ("ID", "experiment", "nSample, "sequence", "resultFile", "resultQoIPath" and "resultInputPath").', 'red', True)
      prt('Exiting', 'red', True)
      sys.exit()

    if self.parent.simulation.objectType == "Simulation_AGATH_0D_IDT":
      self.waitTime = 0.001
    else:
      self.waitTime = 0.01


  def PrepareSamples(self, **kwargs):
    # ----- Initializing the remaining elements -----
    distribution = self.parent.distribution
    if "nSample" in kwargs:
      self.nSample      = kwargs["nSample"]
    if "fixSeed" in kwargs:
      fixSeed           = kwargs["fixSeed"]
    else:
      fixSeed           = False

    # ----- Getting sampling type -----
    temp = self.samplerType.split('_')
    if len(temp) == 2:
      sequence = temp[1]
    else:
      sequence = 'None'

    # ----- Sample cubature first to know its number of samples -----
    # For cubature sampling
    if sequence == 'Clenshaw-Curtis' or sequence == 'SecondFejer':
      dim               = distribution.getDimension()
      xRule             = [sequence] * dim
      self.sampler      = CubatureSampling(self, xRule=xRule, distribution=distribution, xLevel=self.xLevel, isSparse=self.isSparse, verbosity=self.verbosity-1)
      self.sampler.Draw()
      self.nSample      = len(self.sampler.xPoint)
    elif sequence == 'Regular':
      dim               = distribution.getDimension()
      xRule             = [sequence] * dim
      self.sampler      = CubatureSampling(self, xRule=xRule, distribution=distribution, nPointPerDim=self.nSample, verbosity=self.verbosity-1)
      self.sampler.Draw()
      self.nSample      = len(self.sampler.xPoint)
    else:
      #Nothing to be done
      pass

    # ----- Creating or expanding the result file -----
    self.firstSample = HandleResultFile(self.resultFile, self.resultQoIPath, self.resultInputPath, self.nSample, len(self.parent.xVariable), self.isContinued)

    # ----- Sample MC last to know the firstSample if needed -----
    # For cubature sampling
    if sequence == 'Clenshaw-Curtis' or sequence == 'SecondFejer' or sequence == 'Regular':
      pass # Already done for cubature before creating the result file
    # For Monte-Carlo-style sampling
    elif sequence == 'None':
      self.sampler      = MC(self, distribution, self.nSample, self.firstSample, verbosity=self.verbosity-1, fixSeed=fixSeed)
      self.sampler.Draw()
    else:
      self.sampler      = QMC(self, distribution, self.nSample, self.firstSample, sequence=sequence, verbosity=self.verbosity-1)
      self.sampler.Draw()

    # ----- Writing the new sample inputs in the result file -----
    result = h5py.File(self.resultFile, 'r+')
    result[self.resultInputPath][self.firstSample:,:] = self.sampler.xPoint[:,:]
    try:
      xWeight     = self.sampler.xWeight
      hasWeights  = True
    except:
      hasWeights  = False
    if hasWeights:
      try:
        del result[self.resultWeightPath]
      except:
        pass
      result.create_dataset(self.resultWeightPath, data=xWeight)
    result.close()

  #def Run_Old(self):
  #  # ----- Preparing the parallel run -----
  #  counter = 0
  #  xThread = []
  #  xSample = []

  #  # Initializing progbar
  #  if self.verbosity >= 0:
  #    myProgBar = ProgBar('Computing samples...', self.nSample)

  #  for sampleIdx in range(self.firstSample, self.firstSample + self.nSample):
  #    # Creating workdir
  #    workDir = "sample_" + str(sampleIdx)
  #    CopyRunFolder(self.parent.simulation.refFolderPath, workDir)

  #    # Modifying the targets according to the sampled values of the variables
  #    for varIdx in range(len(self.parent.xVariable)):
  #      tempVal = self.sampler.xPoint[sampleIdx - self.firstSample][varIdx]
  #      self.parent.xVariable[varIdx].ModifyTarget(tempVal, workDir)

  #    # Launching computation
  #    if True:  # TODO: Something better for HPC (distributed memory parallelization)? Do I really want to mess with MPI ?
  #      prt("Computing sample " + str(sampleIdx), 'green', self.verbosity - 1)
  #      p = subprocess.Popen(["bash", "RunComputation"], cwd=workDir)
  #      xThread.append(p)
  #      xSample.append(sampleIdx)
  #      counter += 1

  #      while True:
  #        if (self.nSample + self.firstSample - sampleIdx < self.nThread) and (len(xThread) > 0):
  #          shouldContinue = True
  #        else:
  #          shouldContinue = False

  #        for thread in xThread:
  #          index = xThread.index(thread)
  #          poll = thread.poll()

  #          if not poll is None:
  #            # The thread has finished => log the result and make room for another thread
  #            sampleID = xSample[index]
  #            directory = "sample_" + str(sampleID)

  #            QoI = self.parent.simulation.RetrieveQoI(directory)

  #            result = h5py.File(self.resultFile, 'r+')
  #            result[self.resultQoIPath][sampleID] = QoI
  #            result.close()

  #            if not (
  #                    QoI == 0.0):  # TODO: Maybe something more robust here would be useful ? This aimed at keeping the failed computations so that they could be re-run... It is now also being used by "Simulation_To_Keep", perhaps not the best hing to do...
  #              EraseRunFolder(directory)

  #            counter -= 1
  #            del xThread[index]
  #            del xSample[index]
  #            if self.verbosity >= 0:
  #              myProgBar.Update()
  #            shouldContinue = True
  #        if (counter < self.nThread) and not shouldContinue:
  #          break

  #  if self.verbosity >= 0:
  #    myProgBar.Terminate()

  #def ExampleAsynchronousSubprocessesNotOkInThisCaseIThink(self): # TODO: Maybe to test another time ?
  #  from multiprocessing.pool import ThreadPool
  #  import subprocess

  #  def work(sample):
  #    my_tool_subprocess = subprocess.Popen('mytool {}'.format(sample), shell=True, stdout=subprocess.PIPE)
  #    line = True
  #    while line:
  #      myline = my_tool_subprocess.stdout.readline()
  #      # here I parse stdout..

  #  num = None  # set to the number of workers you want (it defaults to the cpu count of your machine)
  #  tp = ThreadPool(num)
  #  for sample in all_samples:
  #    tp.apply_async(work, (sample,))

  #  tp.close()
  #  tp.join()



  def Run(self):
    # ----- Preparing the parallel run -----
    threadCounter = 0
    doneCounter   = 0
    xThread       = []
    xSample       = []
    todoIdx       = range(self.firstSample, self.firstSample + self.nSample)

    # Initializing progbar
    if self.verbosity >= 0:
      myProgBar = ProgBar('Computing samples...', self.nSample)


    # TODO: Write results in batch : Minimize opening/closing time of h5 ?
    #for sampleIdx in range(self.firstSample, self.firstSample + self.nSample):
    while doneCounter < self.nSample: # Loop until all samples have been computed
      """ First part: Check if threads are available and finish work for finished jobs """
      while True: # Wait for a free thread

        #if threadCounter < self.nThread: # Exit right here if we can already launch samples
        #  break

        for thread in xThread:
          index = xThread.index(thread)
          poll = thread.poll()

          if not poll is None:
            # The thread has finished => log the result and make room for another thread
            sampleID = xSample[index]
            directory = "sample_" + str(sampleID)

            QoI = self.parent.simulation.RetrieveQoI(directory)

            result = h5py.File(self.resultFile, 'r+')
            result[self.resultQoIPath][sampleID] = QoI
            result.close()

            if not (QoI == 0.0):  # TODO: Maybe something more robust here would be useful ? This aimed at keeping the failed computations so that they could be re-run... It is now also being used by "Simulation_To_Keep", perhaps not the best hing to do...
              EraseRunFolder(directory)

            threadCounter -= 1
            doneCounter   +=1
            del xThread[index]
            del xSample[index]
            if self.verbosity >= 0:
              myProgBar.Update()

        if threadCounter < self.nThread: # Exit if threads have been freed
          break
        else:
          time.sleep(self.waitTime) # Wait to reduce computational load #TODO: something to calibrate the time of sleep according to computation duration

      """ Second part: Launch the computation of the next sample """
      if len(todoIdx) != 0:
        # Select next sample to compute and remove it from the to do list
        sampleIdx = todoIdx[0]
        del todoIdx[0]

        # Creating workdir
        workDir = "sample_" + str(sampleIdx)
        CopyRunFolder(self.parent.simulation.refFolderPath, workDir)

        # Modifying the targets according to the sampled values of the variables
        for varIdx in range(len(self.parent.xVariable)):
          tempVal = self.sampler.xPoint[sampleIdx-self.firstSample][varIdx]
          self.parent.xVariable[varIdx].ModifyTarget(tempVal, workDir)

        # Launching computation
        if True: # TODO: Something better for HPC (distributed memory parallelization)? Do I really want to mess with MPI ?
          prt("Computing sample " + str(sampleIdx), 'green', self.verbosity-1)
          p = subprocess.Popen(["bash", "RunComputation"], cwd=workDir)
          xThread.append(p)
          xSample.append(sampleIdx)
          threadCounter += 1

    if self.verbosity >= 0:
      myProgBar.Terminate()


  def _InitWithBlock(self, BLOCK):
    keys = BLOCK.GetKeysForBlock()
    self.ID               = BLOCK.GetCharacterFromKey('ID')

    if "sequence" in keys:
      sequence            = BLOCK.GetCharacterFromKey('sequence')
      if sequence == 'None':
        self.samplerType      = 'MC'
        self.nSample          = BLOCK.GetIntegerFromKey('nSample')
        self.resultWeightPath = ''
      elif sequence == 'Clenshaw-Curtis' or sequence == 'SecondFejer':
        self.samplerType      = 'Det_' + sequence
        self.xLevel           = BLOCK.GetIntegerFromKey('xLevel') # TODO: Ensure that this is a list
        self.isSparse         = BLOCK.GetLogicalFromKey('isSparse')
        self.resultWeightPath = BLOCK.GetCharacterFromKey('WeightPath')
      elif sequence == 'Regular':
        self.samplerType      = 'Det_' + sequence
        self.nSample          = BLOCK.GetIntegerFromKey('nSample')
        self.resultWeightPath = BLOCK.GetCharacterFromKey('WeightPath')
      else:
        self.samplerType      = 'QMC_' + sequence
        self.nSample          = BLOCK.GetIntegerFromKey('nSample')
        self.resultWeightPath = ''
    else:
      self.samplerType        = 'MC'
      self.nSample            = BLOCK.GetIntegerFromKey('nSample')
      self.resultWeightPath   = ''

    self.resultFile           = BLOCK.GetCharacterFromKey('resultFile')
    self.resultQoIPath        = BLOCK.GetCharacterFromKey('QoIPath')
    self.resultInputPath      = BLOCK.GetCharacterFromKey('InputPath')
    self.isContinued          = BLOCK.GetLogicalFromKey('isContinued')
    if "outputLevel" in keys:
      self.verbosity          = BLOCK.GetIntegerFromKey('outputLevel')
    else:
      self.verbosity          = -1
    if "nThread" in keys:
      self.nThread            = BLOCK.GetIntegerFromKey('nThread')
    else:
      self.nThread            = 1
