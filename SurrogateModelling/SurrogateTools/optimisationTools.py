import numpy as np
from scipy.optimize import minimize, differential_evolution, dual_annealing, shgo
from scale_rand import *
from Numerics.Optimisation.Dichotomy  import DichotomyOptim     as DOpti

def optimisationTools(fun,strategy, AA=None, b=None, Aeq=None, beq=None, lb=None, ub=None, nonlcon=None, initialGuess=None):
  # Details: Utility function that connects to different optimization schemes
  #
  # inputs:
  # fun - Function handler for optimizable function
  # strategy - Optimization strategy
  # AA - Linear inequality matrix
  # b - Right hand side vector of linear inequality
  # Aeq - Linear equality matrix
  # beq - Linear equality right hand side
  # lb - Lower optimization bound
  # ub - Uppper optimization bound
  # nonlcon - Nonlinear constraint definition
  #
  # outputs:
  # x_opti - Optimized point
  
  #addpath('help_functions')
  
  n = np.size(lb)

  bounds = []
  if lb is not None:
    if ub is not None:
      for i in range(n):
        bounds.append((lb[i], ub[i]))
    else:
      for i in range(n):
        bounds.append((lb[i], 1e30))
  else:
    if ub is not None:
      for i in range(n):
        bounds.append((-1e30, ub[i]))
    else:
      bounds = None

  if strategy == 'evolution' or strategy == 'default':
    x_opti  = differential_evolution(fun, bounds, polish=True, popsize=40, workers=1)
  elif strategy == 'dual_annealing':
    x_opti  = dual_annealing(fun, bounds)
  elif strategy == 'shgo':
    x_opti  = shgo(fun, bounds)
  elif strategy == 'quick':
    #x0      = scale_rand(lb,ub)
    if initialGuess is not None:
      x0 = initialGuess
    else:
      x0 = [3.0, 0.7, 0.7]
    x_opti  = minimize(fun, x0, method=None, bounds=bounds)
  elif strategy == 'multistart':
    xX = []
    xY = []
    for i in range(50):
      x0      = scale_rand(lb,ub)
      temp  = minimize(fun, x0, method=None, bounds=bounds)
      xX.append(temp['x'])
      xY.append(temp['fun'])
    i_opti = np.array(xY).argmin()
    x_opti = xX[i_opti]

  elif strategy == 'GeomDichotomy':
    if len(ub) == 1:
      x_opti  = [DOpti(fun, lb, ub, 10, 1e-3, dichoType='geometric')]
    else:
      print("Dichotomy algorithm unavailable in dimension higher than 1. Reverting to default algorithm with initiation at the lower bound")
      x_opti = minimize(fun, lb, method=None, bounds=bounds)
  else:
    raise Exception('Unknown optimisation strategy : "' + strategy + '"')
  return x_opti
