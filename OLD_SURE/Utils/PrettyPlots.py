import matplotlib.pyplot                as plt
import matplotlib as mpl
font = {'family' : 'sans',
    'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
import math
import numpy as np
import openturns as ot
from Utils.I_O import FunctionSilencer2
from mpl_toolkits.axes_grid1 import AxesGrid

simple  = (7,5)
doubleh = (14, 5)
doublev = (7, 10)
quad    = (14,10)



# TODO: Add functions for summary plot, model plot for 1D and 2D functions, pdf plots...

def AdjustSubplotSpacings(doubleVert=False, doubleHor=False):
  if doubleVert:
    h = 0.7
  else:
    h = 0.35
  if doubleHor:
    w = 0.3
  else:
    w = 0.2
  plt.subplots_adjust(wspace=w, hspace=h)

def CorrectAxSize(ax, fs=20):
  for item in [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels():
    if item is not None:
      item.set_fontsize(fs)

def StripAx(ax, xMin=None, xMax=None, xSpecialX=None, xSpecialY=None):
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  if (xMin is not None) and (xMax is not None):
    if xSpecialX is not None:
      xSpecialX.sort()
      xticks = [xMin[0]] + xSpecialX + [xMax[0]]
    else:
      xticks = [xMin[0]] + [xMax[0]]
    if xSpecialY is not None:
      xSpecialY.sort()
      yticks = [xMin[1]] + xSpecialY + [xMax[1]]
    else:
      yticks = [xMin[1]] + [xMax[1]]
    ax.set_xlim([xMin[0], xMax[0]])
    ax.set_ylim([xMin[1], xMax[1]])
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)


def _normPDF(x):
  return (1.0/(math.sqrt(2.0*math.pi))) * math.exp(-0.5*x*x)

class LegacyKDE:
  def __init__(self, data):
    self.bandwidth = None
    self.data = data.reshape(-1, 1)

  def Evaluate(self,x):
    h    = self.bandwidth
    n = len(self.data)
    temp = 0
    for i in range(n):
        temp += _normPDF((x-self.data[i])/h)
    temp /= h
    temp /= n
    return temp

  def OptimizeBandwidth(self):
    from sklearn.neighbors import KernelDensity
    from sklearn.model_selection import GridSearchCV
    # use grid search cross-validation to optimize the bandwidth
    params = {'bandwidth': np.logspace(-3, 3, 20)}
    grid = GridSearchCV(KernelDensity(), params, cv=5, iid=False)
    grid.fit(self.data)
    self.bandwidth = grid.best_estimator_.bandwidth

  def GetDataPlot(self, nPoint=1000):
    x = np.linspace(self.data.min(), self.data.max(), nPoint)
    y = np.zeros(nPoint)
    for i in range(nPoint):
      y[i] = self.Evaluate(x[i])
    return x, y

  def GetBandwidth(self):
    return self.bandwidth

  def SetBandwidth(self, h):
    self.bandwidth = h

class KDE:
  def __init__(self, data, verbosity=0):
    self.bandwidth=None
    self.result = None
    dataSize = len(data)
    self.data = ot.Sample(dataSize, 1)
    for i in range(dataSize):
      self.data[i,0] = data[i]
    self.kernel = ot.KernelSmoothing()
    self.minData = data.min()
    self.maxData = data.max()
    self.verbosity = verbosity

  def OptimizeBandwidth(self, verySmooth=False):
    if self.verbosity > 0:
      self._InsideOptimizeBandwidth(verySmooth)
    else:
      with FunctionSilencer2():
        self._InsideOptimizeBandwidth(verySmooth)

  def _InsideOptimizeBandwidth(self, verySmooth):
    self.result = self.kernel.build(self.data)
    self.bandwidth = self.kernel.getBandwidth()[0]
    if verySmooth:
      self.result = self.kernel.buildAsKernelMixture(self.data, ot.Point([2.0 * self.bandwidth]))
  #def OptimizeBandwidth(self, verySmooth):
  #  self.result = self.kernel.build(self.data)
  #  self.bandwidth = self.kernel.getBandwidth()[0]
  #  if verySmooth:
  #    self.result = self.kernel.buildAsKernelMixture(self.data, ot.Point([2.0 * self.bandwidth]))

  def Evaluate(self,x):
    return self.result.computePDF(x)

  def GetDataPlot(self, nPoint=1000):
    x = np.linspace(self.minData, self.maxData, nPoint)
    y = np.zeros(nPoint)
    for i in range(nPoint):
      y[i] = self.Evaluate(x[i])
    return x,y

  def GetBandwidth(self):
    return self.bandwidth




def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.
    credits to Paul H: https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
