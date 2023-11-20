import numpy as np
from Utils.I_O import ErrorMsg, WarningMsg
from Sampling.Deterministic import QuadratureSampling
from Sampling.Deterministic import CubatureSampling
import sys

fileName = 'Numerics/Integration.py'

def Integrate(func, bounds, **kwargs):
  # func is a link to a python function that takes a vector or scalar as input and outputs a scalar. Can be an analytical function or something that reads a csv if the quadrature points are already known
  # bounds is a list of two numbers for 1D integration or a list of lists of two scalars for n-dimensional integration
  
  if isinstance(bounds[0], list):
    return IntegrateND(func, bounds, kwargs)
  else:
    return Integrate1D(func, bounds, kwargs)

def Integrate1D(func, bounds, kwargs):
  # Prepare quadrature points according to the inputs
  if 'method' in kwargs:
    method = kwargs['method']
  else:
    ErrorMsg('A "method" argument must be provided', fileName)
  if method in ['Rectangle', 'Trapeze']:
    if 'nPoint' in kwargs:
      nPoint = kwargs['nPoint']
    else:
      ErrorMsg('For "Rectangle" and "Trapeze" methods, the argument "nPoint" must be provided', fileName)
    Drawer = QuadratureSampling(None, 'Regular', 'Uniform', bounds, nPoint=nPoint)
  elif method in ['SecondFejer', 'Clenshaw-Curtis', 'Simpson']:
    if 'level' in kwargs:
      level = kwargs['level']
    else:
      ErrorMsg('For "SecondFejer", "Clenshaw-Curtis" and Simpson methods, the argument "level" must be provided', fileName)
    Drawer = QuadratureSampling(None, method, 'Uniform', bounds, level=level)
  else:
      ErrorMsg('Invalid "method" argument. Valid arguments are "Rectangle", "Trapeze", "SecondFejer", "Clenshaw-Curtis" and "Simpson"', fileName)

  # Computing the integral
  Drawer.Draw()
  result = 0.0
  if method in ['Rectangle', 'SecondFejer', 'Clenshaw-Curtis', 'Simpson']:
    for i in range(len(Drawer.xPoint)):
      result += func(Drawer.xPoint[i]) * Drawer.xWeight[i]
  else: # Trapeze method
    eval1 = func(Drawer.xPoint[0]) * Drawer.xWeight[0] * len(Drawer.xPoint)
    for i in range(1,len(Drawer.xPoint)):
      eval2 = func(Drawer.xPoint[i]) * Drawer.xWeight[i] * len(Drawer.xPoint)
      result += (eval2 + eval1) / 2
      eval1 = eval2
    result /=(len(Drawer.xPoint)-1)

  return result


def IntegrateND(func, bounds, kwargs):
  # Prepare quadrature points according to the inputs
  if 'samplingMethod' in kwargs:
    samplingMethod = kwargs['samplingMethod']
  else:
    ErrorMsg('A "samplingMethod" argument must be provided', fileName)

  if samplingMethod == 'Cubature':
    if 'isSparse' in kwargs:
      isSparse = kwargs['isSparse']
    else:
      isSparse = False
    if 'quadratureMethod' in kwargs:
      quadMethod = kwargs['quadratureMethod']
    else:
      ErrorMsg('For "Cubature" samplingMethod, a "quadratureMethod" argument must be provided', fileName)
    if quadMethod in ['SecondFejer', 'Clenshaw-Curtis', 'Simpson']:
      if 'maxLevel' in kwargs:
        maxLevel = kwargs['maxLevel']
      else:
        ErrorMsg('For "SecondFejer", "Clenshaw-Curtis" and "Simson" methods, the argument "maxLevel" must be provided', fileName)
      Drawer = CubatureSampling(None, [quadMethod]*len(bounds), ['Uniform']*len(bounds), bounds, maxLevel=maxLevel, isSparse=isSparse)
      Drawer.Draw()
    else:
        ErrorMsg('Invalid "quadratureMethod" argument. Valid arguments are "SecondFejer", "Clenshaw-Curtis" and "Simpson"', fileName)
  elif samplingMethod == 'MC':
    ErrorMsg('"MC" integration has not yet been implemented. Please come back later', fileName)
  elif samplingMethod == 'QMC':
    ErrorMsg('"QMC" integration has not yet been implemented. Please come back later', fileName)
  else:
    ErrorMsg('Invalid "samplingMethod" argument. Valid arguments are "Cubature", "MC" and "QMC"', fileName)

  # Computing the integral
  result = 0.0
  for i in range(len(Drawer.xPoint)):
    result += func(Drawer.xPoint[i]) * Drawer.xWeight[i]
  if 'plot' in kwargs:
    display = kwargs['plot']
  else:
    display = False
  if display:
    if len(bounds)==2:
      import matplotlib.pyplot as plt
      x = np.zeros(np.size(Drawer.xPoint, 0))
      x[:] = Drawer.xPoint[:,0]
      y = np.zeros(np.size(Drawer.xPoint, 0))
      y[:] = Drawer.xPoint[:,1]
      fig = plt.figure(figsize=(14,10))
      plt.scatter(x, y, label='Cubature')
      plt.legend(loc='best')
      plt.show()
    else:
      WarningMsg('The plot of cubature points is not available for dimensions superior than 2', fileName, True)

  return result
