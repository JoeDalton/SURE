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

xTrue2D = []
with open('DIST2D_77.tidat', 'r') as f:
  lines = f.readlines()
  for line in lines:
    xTrue2D.append(float(line))
xTrue2D = np.array(xTrue2D)[:]

##Dataset options
#step = 1/64.0
#level = 2
#dim = 2
#xMin = -2.0 * step
#xMax = 2.0 * step

## Build dummy function for tests
#intRule = 'Simpson'
#distributionType = 'Uniform'
#xDistributionParam = [-1.0, 1.0]
#dim=2
#level=2
#drawer = CubatureSampling(None, xRule=[intRule]*dim, xDistributionType=[distributionType]*dim, xxDistributionParam = [xDistributionParam]*dim, xLevel = [level]*dim, verbosity = 0)
#drawer.Draw()
#xSample = drawer.xPoint

#def dummyFunction(point):
#  return point[0]*point[0]+point[1]*point[1]
#for i in range(len(xSample)):
#  xTrue2D[i] = dummyFunction(xSample[i])

level=3
dim = 2
intRule = 'Simpson'
distributionType = 'Uniform'
xDistributionParam = [0.0, 1.0]
idrawer = CubatureSampling(None, xRule=[intRule]*dim, xDistributionType=[distributionType]*dim, xxDistributionParam = [xDistributionParam]*dim, xLevel = [level]*dim, verbosity = 0)
idrawer.Draw()
xiSample = idrawer.xPoint



from Legendre import LegendreUnivariate
from Projection import *
myCustomPCE = BuildPCE(2,3, xTrue2D)


def customPCE(x, xMin, xMax):
  #Entre -1 et 1 !!!
  xNorm2 = 2.0 * (x - xMin) / (xMax - xMin) - 1.0
  return myCustomPCE(xNorm2)




from OK_model import *
myCustomKrig = OK_model("Matern32", xiSample, xTrue2D, opti="default")
#
##for i in range(nReSample*10):
##  xCustomKrigEval[i] = myCustomKrig.predict(xX[i])[0]















#factor = 2.0/(4.0*step)

#print(myCustomPCE.EvaluateCurvature(0)*factor*factor)

#dummy:
#print(myCustomPCE.EvaluateCurvature([0.0,0.0]))
#quit()
##
##
##print(myCustomPCE.EvaluateGradient([0.0,0.0]))
##xTrue2D = np.reshape(xTrue2D, ((5,5)))
##fxFD = (xTrue2D[2,3] - xTrue2D[2,1])/1.0
##fyFD = (xTrue2D[3,2] - xTrue2D[1,2])/1.0
##print(fxFD, fyFD)
#quit()


"""
Plot stuff
"""
##Surrogate options
dim =2
intRule = 'Simpson'
distributionType = 'Uniform'
xDistributionParam = [0.0, 1.0]
drawer = CubatureSampling(None, xRule=[intRule]*dim, xDistributionType=[distributionType]*dim, xxDistributionParam = [xDistributionParam]*dim, xLevel = [level]*dim, verbosity = 0)
drawer.Draw()
xSample = drawer.xPoint
#temp = np.reshape(xSample, ((5,5,2)))
#xSample = np.reshape(temp, ((25,2)), order='F')
temp = np.reshape(xSample, ((7,7,2)))
xSample = np.reshape(temp, ((49,2)), order='F')
#xWeight = drawer.xWeight
#xEval = xTrue2D
#
#
## Building PCE
#objPCE = PCE(ID='titi', verbosity=1, distribution=ot.ComposedDistribution([ot.Uniform(0.0,1.0)]*dim))
#objPCE.verbosity = True
#objPCE.DefinePolyBasis() 
#objPCE.DefineTruncatureStrat(strategy="Fixed", maxTotDegree = polDegree)
#objPCE.DefineEvalStrat(evalStrategy="Cubature")
#objPCE.Build(xSample, xEval, xWeight=xWeight)
#print(objPCE.polyCoeffs)
#print(objPCE.responseSurface)
#
#def myPCE(x, xMin, xMax):
#  xNorm = np.zeros_like(x)
#  xNorm[:] = (x[:] - xMin) / (xMax - xMin)
#  return objPCE.Evaluate([xNorm])[0]


#############################################################
# STEP 4: Some plots
#############################################################
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
# DEBUG. TODO: better thingy here
xMin = -1
xMax = 1
yMin = xMin
yMax = xMax
nPt = 50
x = np.linspace(xMin, xMax, nPt)
y = np.linspace(yMin, yMax, nPt)
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)
ZK = np.zeros_like(X)
for i in range(nPt):
  for j in range(nPt):
    #Z[i, j] = myPCE(np.array([X[i, j], Y[i, j]]), xMin, xMax)[0]
    Z[i, j] = customPCE(np.array([X[i, j], Y[i, j]]), xMin, xMax)
    ZK[i, j] = myCustomKrig.predict(np.array([X[i, j], Y[i, j]]))[0]

norm = plt.Normalize(Z.min(), Z.max())
colors = cm.hot(norm(Z))
rcount, ccount, _ = colors.shape

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, Z, cmap=cm.coolwarm, alpha=0.3, label='True')
ax.scatter(xMin + xSample[:,0] * (xMax-xMin), xMin + xSample[:,1] * (xMax-xMin), xTrue2D[:], marker='o', color='r', label='Reference samples')
#ax.plot_wireframe(xMin + X * (xMax-xMin), xMin + Y * (xMax-xMin), ZK,color='k', label='Krig')
ax.plot_wireframe(yMin + Y * (yMax-yMin), xMin + X * (xMax-xMin), ZK,color='k', label='Kriging')
#ax.plot_wireframe(X,Y, ZK,color='k', label='Krig')
ax.view_init(azim=0, elev=90)
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
pp.CorrectAxSize(ax)
plt.legend(loc='best')
#plt.savefig(modFile, bbox_inches='tight') #DEBUG
plt.show()
