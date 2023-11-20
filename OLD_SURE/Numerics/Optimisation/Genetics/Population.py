import random
random.seed(0)
import numpy    as np
import multiprocessing
from math       import floor #*
from Numerics.Optimisation.Genetics.Eval  import MultiEval
from operator   import attrgetter
from tqdm       import tqdm
import sys
from Utils.I_O import ErrorMsg


########################################################
########################################################
fileName = 'Numerics/Optimisation/Genetics/Population.py'

##############################
##          LAB RAT         ##
##############################


class LabRat():
  # Properies
  xGen        = []
  parent      = None
  grade       = None

  # Methods
  def __init__(self, parent):
    self.grade        = None
    self.parent       = parent

  def Evaluate(self):
    self.grade        = self.parent.GradingFunc(self.xGen)
  
  def EvaluateForMulti(self):
    self.grade        = self.parent.GradingFunc(self.xGen)
    return self.grade

  def BreedRandom(self):
    xLBound = self.parent.xLBound
    xUBound = self.parent.xUBound
    temp = np.random.uniform(0.0,1.0,len(xLBound))
    self.xGen = xLBound + temp * (xUBound - xLBound)
 
  def BreedInterpolate(self, mIdx, fIdx, weight=0.5):
    mXGen = self.parent.xRat[mIdx].xGen
    fXGen = self.parent.xRat[fIdx].xGen
    self.xGen = weight*mXGen + (1-weight)*fXGen
 
  def BreedCrossover(self, mIdx, fIdx):
    mXGen = self.parent.xRat[mIdx].xGen
    fXGen = self.parent.xRat[fIdx].xGen
    self.xGen = np.zeros_like(mXGen)
    for i in range(len(self.xGen)):
      if random.random()>0.5:
        self.xGen[i] = mXGen[i]
      else:
        self.xGen[i] = fXGen[i]
 
  def BreedMutate(self, mIdx, scale):
    mXGen     = self.parent.xRat[mIdx].xGen
    #rand      = np.ones_like(mXGen) + np.random.uniform(-scale,scale,len(mXGen))
    realScale = scale * (self.parent.xUBound - self.parent.xLBound)
    tempRand  = np.ones_like(mXGen) + np.random.normal(loc=0.0, scale=1.0, size=len(mXGen))
    rand      = realScale * tempRand
    self.xGen = mXGen + rand
    self.CheckRange()
  
  def CheckRange(self):
    self.xGen = np.clip(self.xGen, self.parent.xLBound, self.parent.xUBound)

  #def proximity(gen1,gen2):
  #    proximity = 0
  #    nGen = len(gen1)
  #    for i in range(nGen):
  #        d = gen1[i] - gen2[i]
  #        proximity+= d*d
  #    return sqrt(proximity)

  #def gen2log(gen,logfile):
  #  #logfile.write(" %3.2f          " %gen[0])
  #  logfile.write (' {0:10.2f}       '.format(gen[0]))
  #  nGen = len(gen[1])
  #  for i in range(nGen):
  #    logfile.write(" %1.4f " %gen[1][i])
  #  logfile.write("\n")




  
########################################################
########################################################


##############################
##        POPULATION        ##
##############################


class Population():
  # Properties
  nRat            = 0
  xRat            = []
  xLBound         = None # np array
  xUBound         = None # np array
  GradingFunc     = None # function
  keepFrac        = 0.0
  interProba      = 0.0
  crossProba      = 0.0
  mutatProba      = 0.0
  no1ParentProba  = 0.0
  mutationScale   = 0.0
  isFirstEval     = True

  # Methods
  def __init__(self, xLBound, xUBound, GradingFunc, **kwargs): #keepFrac, interProba, crossProba, mutatProba, no1ParentProba, mutationScale):
    self.xRat           = []
    self.GradingFunc    = GradingFunc
    self.xLBound        = np.array(xLBound)
    self.xUBound        = np.array(xUBound)
    self.isFirstEval    = True
    if "preSet" in kwargs:
      preSet = kwargs["preSet"]
      if preSet == "I'm Too Young To Die":
        self.nRat           = 100
        self.keepFrac       = 0.3
        self.interProba     = 0.2
        self.crossProba     = 0.1
        self.mutatProba     = 0.2
        self.no1ParentProba = 0.5
        self.mutationScale  = 0.1
      elif preSet == "Hurt Me Plenty":
        self.nRat = 1000
        self.keepFrac = 0.4
        self.interProba = 0.2
        self.crossProba = 0.1
        self.mutatProba = 0.2
        self.no1ParentProba = 0.4
        self.mutationScale = 0.1
      elif preSet == "Nightmare":
        self.nRat = 10000
        self.keepFrac = 0.5
        self.interProba = 0.2
        self.crossProba = 0.1
        self.mutatProba = 0.2
        self.no1ParentProba = 0.2
        self.mutationScale = 0.1
      elif preSet == "Hold The Line":
        self.nRat = 1000
        self.keepFrac = 0.1
        self.interProba = 0.01
        self.crossProba = 0.0
        self.mutatProba = 0.8
        self.no1ParentProba = 0.8
        self.mutationScale = 0.05
      else:
        ErrorMsg('Unsupported preset: "' + preSet + '"', fileName)
    else:
      print("TODO: custom init")
      sys.exit()
    #TODO: Check that inter + cross + mutat < 1



  #def initBkp(self, nRat, xLBound, xUBound, GradingFunc, keepFrac, interProba, crossProba, mutatProba, no1ParentProba, mutationScale):
  #  self.nRat           = nRat
  #  self.xRat           = []
  #  self.xLBound        = np.array(xLBound)
  #  self.xUBound        = np.array(xUBound)
  #  self.GradingFunc    = GradingFunc
  #  self.keepFrac       = keepFrac
  #  self.interProba     = interProba
  #  self.crossProba     = crossProba
  #  self.mutatProba     = mutatProba
  #  self.no1ParentProba = no1ParentProba
  #  self.mutationScale  = mutationScale
  #  self.isFirstEval    = True
  #  #TODO: Check that inter + cross + mutat < 1

  def Populate(self, labRats=[]):
    if len(labRats)<self.nRat:
      self.xRat = labRats
    else:
      self.xRat = labRats[:self.nRat]
    for i in range(self.nRat-len(labRats)):
      myRat = LabRat(self)
      myRat.BreedRandom()
      self.xRat.append(myRat)
  
  def Evaluate(self):
    if self.isFirstEval:
      nBest = 0
      self.isFirstEval = False
    else:
      nBest = int(floor(self.nRat * self.keepFrac)) # We only evaluate the new rats
    for i in range(nBest,self.nRat):
      self.xRat[i].Evaluate()
    self.Sort()

  def EvaluateParallel(self):
    if self.isFirstEval:
      nBest = 0
      self.isFirstEval = False
    else:
      nBest = floor(self.nRat * self.keepFrac) # We only evaluate the new rats
    p = multiprocessing.Pool(multiprocessing.cpu_count())
    result = p.map(self.SingleEval, tqdm(range(nBest,self.nRat)))
    p.close()
    p.join()

    for i in range(nBest,self.nRat):
      self.xRat[i].grade = result[i-nBest] #TODO: vectorized expression to speed up ?
    self.Sort()

  def SingleEval(self,i):
    return self.xRat[i].EvaluateForMulti()

  def Sort(self):
    self.xRat = sorted(self.xRat, key=attrgetter('grade'))
    #self.xRat[0].grade.Print()

  def GetBestGrade(self):
    # We assume the population is already sorted (end of eval_fitness)
    return self.xRat[0].grade #.xGrade[-1]

  def Evolve(self):
    newXRat = []
    nBest = int(floor(self.nRat * self.keepFrac))
    # Keep the best rats for the next gen
    for i in range(nBest):
      newXRat.append(self.xRat[i])
    
    # Replace the others with offsprings
    for i in range(nBest,self.nRat): # TODO: Parallelize this loop ? Not critical ?
      # Initialize new child
      myRat = LabRat(self)

      # Choose parents
      var1 = random.random()
      if var1 < self.no1ParentProba:
        mIdx = 0                          # Choose the two best rats to reproduce
        fIdx = 1
      else: 
        mIdx = random.randint(0,nBest-1)  # Choose any 2 in the nBest
        fIdx = random.randint(0,nBest-1)
      
      # Choose reproduction method
      var2 = random.random()
      if var2 < self.interProba:                                        # Interpolation
        myRat.BreedInterpolate(mIdx, fIdx)
      elif var2 < self.interProba + self.crossProba:                    # Crossover
        myRat.BreedCrossover(mIdx, fIdx)
      elif var2 < self.interProba + self.crossProba + self.mutatProba:  # Mutation
        myRat.BreedMutate(mIdx, self.mutationScale)
      else:                                                             # Random
        myRat.BreedRandom()

      # Add the new child to the population
      newXRat.append(myRat)

    # Replace the old population by the new
    self.xRat = newXRat
