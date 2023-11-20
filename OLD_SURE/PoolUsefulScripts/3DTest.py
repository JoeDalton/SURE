#!/usr/bin/env python2.7

from __future__ import print_function
import h5py
import time
import Link
import numpy                                      as np
from numpy import sin
import openturns                                  as ot
import matplotlib.pyplot                          as plt
from SurrogateModelling.PCE             import *
import Utils.PrettyPlots                          as pp
from Sampling.Deterministic import *

t0 = time.time()

xTrue3D = []
with open('DIST3D_sample1.tidat', 'r') as f:
  lines = f.readlines()
  for line in lines:
    xTrue3D.append(float(line))
xTrue3D = np.array(xTrue3D)[:]

from Legendre import LegendreUnivariate
from Projection import *
myCustomPCE = BuildPCE(3, xTrue3D)


xTrue3D = np.reshape(xTrue3D, ((5,5,5)))
print(xTrue3D)
def EvaluateCurvature():
  dim = 3
  one = 1.0
  two = 2.0
  epsilon = 1e-15
  fx    = xTrue3D[3,2,2] - xTrue3D[1,2,2]
  fy    = xTrue3D[2,3,2] - xTrue3D[2,1,2]
  if dim == 3:
    fz  = xTrue3D[2,2,3] - xTrue3D[2,2,1]
  fxx   = (xTrue3D[4,2,2] - xTrue3D[2,2,2]) - (xTrue3D[2,2,2] - xTrue3D[0,2,2])
  fyy   = (xTrue3D[2,4,2] - xTrue3D[2,2,2]) - (xTrue3D[2,2,2] - xTrue3D[2,0,2])
  fxy   = (xTrue3D[3,3,2] - xTrue3D[3,1,2]) - (xTrue3D[1,3,2] - xTrue3D[1,1,2])
  if dim == 3:
    fzz = (xTrue3D[2,2,4] - xTrue3D[2,2,2]) - (xTrue3D[2,2,2] - xTrue3D[2,2,0])
    fzx = (xTrue3D[3,2,3] - xTrue3D[3,2,1]) - (xTrue3D[1,2,3] - xTrue3D[1,2,1])
    fyz = (xTrue3D[2,3,3] - xTrue3D[2,3,1]) - (xTrue3D[2,1,3] - xTrue3D[2,1,1])

  if dim == 2:
    norm  = math.sqrt( fx**2 + fy**2)
    nx    = fx/norm
    ny    = fy/norm
    kappa = -((one-nx**2)*fxx+(one-ny**2)*fyy- \
                                 two*nx*ny*fxy)/max(norm,epsilon)
  elif dim ==3:
    norm  = math.sqrt( fx**2 + fy**2 + fz**2 )
    nx    = fx/norm
    ny    = fy/norm
    nz    = fz/norm
    kappa = -((one-nx**2)*fxx+(one-ny**2)*fyy+(one-nz**2)*fzz- \
                                 two*nx*ny*fxy-two*nx*nz*fzx-two*ny*nz*fyz)/ \
                                 max(norm,epsilon)
  else:
    raise Exception("Curvature error, this should not happen")
  return kappa

print(myCustomPCE.EvaluateCurvature([0.0,0.0, 0.0]))
print(EvaluateCurvature())
