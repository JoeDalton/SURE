#!/usr/bin/env python2.7

from __future__ import print_function
import h5py
import time
import Link
import numpy                                      as np
from numpy import sin,cos
import openturns                                  as ot
import matplotlib.pyplot                          as plt
from SurrogateModelling.PCE             import *
import Utils.PrettyPlots                          as pp
from Sampling.Deterministic import *

t0 = time.time()



#Dataset options
level = 2
dim = 1
xMin = -2.0
xMax = 2.0

def myFunc(x, xMin, xMax):
  return cos(x+1)

def myFuncNorm(xNorm, xMin, xMax):
  x = xMin + xNorm * (xMax - xMin)
  return myFunc(x, xMin, xMax)



#Surrogate options
polDegree =2
intRule = 'Simpson'
distributionType = 'Uniform'
xDistributionParam = [0.0, 1.0]

#Plotting options
nReSample = 1000


# Sampling
drawer = CubatureSampling(None, xRule=[intRule], xDistributionType=[distributionType], xxDistributionParam = [xDistributionParam], xLevel = [level], verbosity = 0)
drawer.Draw()
xSample = drawer.xPoint
xWeight = drawer.xWeight
xEval = myFuncNorm(xSample[:,0], xMin, xMax)

print(xSample)
#quit()


# Building PCE
objPCE = PCE(ID='titi', verbosity=1, distribution=ot.ComposedDistribution([ot.Uniform(0.0,1.0)]*dim))
objPCE.verbosity = True
objPCE.DefinePolyBasis() 
objPCE.DefineTruncatureStrat(strategy="Fixed", maxTotDegree = polDegree)
objPCE.DefineEvalStrat(evalStrategy="Cubature")
objPCE.Build(xSample, xEval, xWeight=xWeight)

xX        = np.linspace(xMin*10, xMax*10, nReSample*10)
xPCEEval  = np.zeros(nReSample*10)
xCustomPCEEval  = np.zeros(nReSample*10)
xCustomKrigEval  = np.zeros(nReSample*10)
xCustomKrigLC    = np.zeros(nReSample*10)
xCustomKrigUC    = np.zeros(nReSample*10)
xTrueEval = myFunc(xX, xMin, xMax)




def myPCE(x, xMin, xMax):
  xNorm = (x - xMin) / (xMax - xMin)
  return objPCE.Evaluate([xNorm])[0]

for i in range(nReSample*10):
  #xPCEEval[i] = objPCE.Evaluate([xX[i]])[0]
  xPCEEval[i] = myPCE(xX[i], xMin, xMax)



def slopeFDOT(x, xMin, xMax):
  deltaX = 0.00005
  res = (myPCE(x+deltaX, xMin, xMax) - myPCE(x-deltaX, xMin, xMax))/(2*deltaX)
  return res

#print(slopeFDOT(0,xMin, xMax))


#"""
#PCE by hand
#"""
#from Legendre import LegendreUnivariate
#from Projection import *
#myCustomPCE = BuildPCE(1, xEval)
#
#
#def customPCE(x, xMin, xMax):
#  #Entre -1 et 1 !!!
#  xNorm2 = 2.0 * (x - xMin) / (xMax - xMin) - 1.0
#  return myCustomPCE(xNorm2)
#
#for i in range(nReSample*10):
#  xCustomPCEEval[i] = customPCE(xX[i], xMin, xMax)

"""
Kriging by hand
"""
from OK_model import *
#myCustomKrig = OK_model("Matern32", xSample, xEval, "default", nugget=1.0)
#myCustomKrig = OK_model("Matern32", xSample, xEval, "default", givenTheta=[3.0])
myCustomKrig = OK_model("Matern52", xSample, xEval, opti="default", givenTheta=[3.0])
#myCustomKrig = OK_model("Exp", xSample, xEval, "default", givenTheta=[3.0])
#myCustomKrig = OK_model("SqExp", xSample, xEval, opti="default")
#myCustomKrig = OK_model("CubicSpline", xSample, xEval, "default")
print(myCustomKrig.theta)
for i in range(nReSample*10):
  temp = myCustomKrig.predict(xX[i], conf=95)
  xCustomKrigEval[i] = temp[0]
  xCustomKrigLC[i] = temp[2]
  xCustomKrigUC[i] = temp[3]

#testValue = 1
#print("P(x)")
#print(myPCE(testValue, xMin, xMax))
#print(customPCE(testValue, xMin, xMax))




#quit()
#############################################################
# STEP 4: Some plots
#############################################################
fig = plt.figure(1, figsize=(7,5))
pp.AdjustSubplotSpacings()

ax1 = fig.add_subplot(111)
#ax1.plot(xX, xPCEEval, color='k', label='PCE from OT') 
#ax1.plot(xX, xCustomPCEEval+0.01, color='r', label='custom PCE') 
ax1.plot(xMin + xX*(xMax-xMin), xCustomKrigEval, color='b', label='custom Kriging') 
ax1.plot(xX,xTrueEval, color='g', label='True')
ax1.scatter(xMin + xSample*(xMax-xMin),xEval, color='r', label='Samples')

ax1.fill_between(xMin + xX*(xMax-xMin), xCustomKrigLC, xCustomKrigUC, alpha=0.2)

ax1.legend()
ax1.set_xlim(-10,10)
ax1.set_ylim(-4,4)
plt.show()
