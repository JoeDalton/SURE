from Numerics.Optimisation.Genetics.Population import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from Utils.I_O import ErrorMsg
font = {'family' : 'sans',
    'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)

fileName = 'Numerics/Optimisation/Genetics/Genetics.py'

class GeneticAlgorithm():
  # Properties
  population  = None
  nGeneration = 0
  outLevel    = 0
  
  
  # Methods
  def __init__(self, population, outLevel=0, **kwargs):
    self.population   = population
    self.outLevel     = outLevel
    if "nGen" in kwargs:
      self.nGeneration = kwargs["nGen"]
    elif "preSet" in kwargs:
      preSet = kwargs["preSet"]
      if preSet == "I'm Too Young To Die":
        self.nGeneration = 200
      elif preSet == "Hurt Me Plenty":
        self.nGeneration = 2000
      elif preSet == "Nightmare":
        self.nGeneration = 20000
      elif preSet == "Hold The Line":
        self.nGeneration = 2000
      else:
        ErrorMsg('Unsupported preset: "' + preSet + '"', fileName)

  def Run(self, parallel=False, labRats=[]):
    xIte        = []
    xBestGrade  = []
    #fig, sc     = self.PreparePlot(xIte, xBestGrade)

    self.Prt('')
    self.Prt(' ===== Generation: 0')
    xIte.append(0)
    self.population.Populate(labRats=[])
    if parallel:
      self.population.EvaluateParallel()
    else:
      self.population.Evaluate()
    xBestGrade.append(self.population.GetBestGrade())
    #self.UpdatePlot(fig, sc, xIte, xBestGrade)
    self.Prt("Best Guy: " + str(self.population.xRat[0].xGen))
    self.Prt("Result: " + xBestGrade[-1].Display())

    for i in range(1, self.nGeneration):
      #TODO 2: Add phases of exploration/exploitation: speciation ? Or restart algorithm a few times after convergence
      self.Prt('')
      self.Prt(' ===== Generation: ' + str(i))
      xIte.append(i)
      self.population.Evolve()
      if parallel:
        self.population.EvaluateParallel()
      else:
        self.population.Evaluate()
      xBestGrade.append(self.population.GetBestGrade())
      #self.UpdatePlot(fig, sc, xIte, xBestGrade)
      self.Prt("Best Guy: " + str(self.population.xRat[0].xGen))
      self.Prt("Result: " + xBestGrade[-1].Display())

      # Exit condition:
      if i > 50:
        if xBestGrade[i] >= xBestGrade[i-50]:
          self.Prt('')
          self.Prt('Convergence reached')
          break

  def Prt(self, message):
    if self.outLevel > 0:
      print(message)

  def PreparePlot(self, xIte, xBestGrade):
    fig = plt.figure(figsize=(10,6))
    if self.outLevel > 2:
      plt.ion() # Interactive plot
      ax1 = fig.add_subplot(111)
      xIte=[]
      sc = ax1.scatter(xIte, xBestGrade, color='k', marker='+', label='Best grade')    
      ax1.set_xlabel(r'Iteration')
      ax1.set_ylabel(r'Grade')
      ax1.set_xlim([0.0,nIte])
      ax1.set_ylim([1e-5,1])
      plt.legend(loc='best')
      plt.yscale('log')
      plt.grid(True, axis='y')
      for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels():
        item.set_fontsize(18)
      return fig, sc
    else:
      return None, None

  def UpdatePlot(self, fig, sc, xIte, xBestGrade):
    if self.outLevel > 2:
      sc.set_offsets(np.c_[xIte,xBestGrade])
      fig.canvas.draw_idle()
      plt.pause(0.1)

  def KeepPlot(self, fig, sc):
    if self.outLevel > 2:
      plt.show(block=True)




