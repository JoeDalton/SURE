#!/usr/bin/env python2.7

import Link
import numpy                as np
import matplotlib.pyplot    as plt
import openturns            as ot
from Sampling.MC            import *
from Sampling.QMC           import *
from Sampling.RQMC          import *
from Sampling.Deterministic import *
import Utils.PrettyPlots as pp


#------------------------------------------------------------ Preparation of samples
#nSample      = 1000
nSample      = 10
firstSample  = 0
nQMCDrawer   = 10
#distribution = ot.ComposedDistribution([ot.Uniform(0.0, 1.0)] * 2)
distribution = ot.ComposedDistribution([ot.Uniform(0.0, 1.0),ot.Normal(0.0, 1.0)])
verb = -1


############### MC
MCDrawer   = MC(None, distribution, nSample, firstSample=firstSample, verbosity=verb, fixSeed=False)
s1 = MCDrawer.Draw()
x1 = np.zeros(nSample)
y1 = np.zeros(nSample)
x1 = s1[:, 0]
y1 = s1[:, 1]

############### QMC
QMCDrawer  = QMC(None, distribution, nSample, firstSample=firstSample, sequence='Sobol', verbosity=verb)
s2 = QMCDrawer.Draw()
y2 = np.zeros(nSample)
x2 = np.zeros(nSample)
x2 = s2[:, 0]
y2 = s2[:, 1]

############### RQMC
RQMCDrawer = RQMC(None, distribution, nSample/nQMCDrawer, nQMCDrawer=nQMCDrawer, firstSample=firstSample, sequence='Sobol', verbosity=verb, fixSeed=True)
s3 = RQMCDrawer.Draw()
x3 = np.zeros(nSample)
y3 = np.zeros(nSample)
for i in range(len(s3)):
  for j in range(nSample/nQMCDrawer):
    x3[(i*nSample/nQMCDrawer)+j] = s3[i][j, 0]
    y3[(i*nSample/nQMCDrawer)+j] = s3[i][j, 1]

################ Quadrature
#---Nested
#rule               = 'Clenshaw-Curtis'
rule               = 'SecondFejer'
#rule               = 'Simpson'
maxLevel           = 3
#---Regular
#rule               = 'Regular'
#nPoint             = 8

distributionType   = 'Uniform'
xDistributionParam = [0.0, 1.0]
xLevel = range(0, maxLevel)
xQuad = []
x4    = []
y4    = []
for level in xLevel:
  #quadDrawer = QuadratureSampling(None, rule, distributionType, xDistributionParam, level=level, nPoint = nPoint, verbosity=verb)
  #quadDrawer.Draw()
  #x4Temp = quadDrawer.xPoint.tolist()
  #y4Temp = quadDrawer.xPoint.tolist()
  CubDrawer = CubatureSampling(None, xRule=[rule], xDistributionType=[distributionType], xxDistributionParam=[xDistributionParam], xLevel=[level], isSparse=True, verbosity = verb+1) #DEBUG TEST for 1D cubature (= quadrature : yes ! Something fishy for Simpson ?)
  CubDrawer.Draw()
  x4Temp = CubDrawer.xPoint.tolist()
  y4Temp = CubDrawer.xPoint.tolist()
  for k in range(len(y4Temp)):
    y4Temp[k] = level
  x4 = x4 + x4Temp
  y4 = y4 + y4Temp

################## Cubature
dim                 = 2
#rule                = 'Clenshaw-Curtis'
rule                = 'SecondFejer'
distributionType    = 'Uniform'
xDistributionParam  = [0.0,1.0]
level               = 2
xRule               = [rule] * dim
xDistributionType   = [distributionType] * dim
xxDistributionParam = [xDistributionParam] * dim
#xLevel              = [level] * dim
xLevel              = [level, level]

# Full
#CubDrawer = CubatureSampling(None, xRule, xDistributionType, xxDistributionParam, xLevel=xLevel, isSparse = False, verbosity = True)
#CubDrawer = CubatureSampling(None, xRule, xDistributionType, xxDistributionParam, maxLevel=level, isSparse = False, verbosity = True)

#Smolyak
#CubDrawer = CubatureSampling(None, xRule=xRule, xDistributionType=xDistributionType, xxDistributionParam=xxDistributionParam, xLevel=xLevel, isSparse = False, verbosity = verb+1)
CubDrawer = CubatureSampling(None, rule=rule, distribution=distribution, level=level, isSparse = True, verbosity = verb+1)
CubDrawer.Draw()

x5 = np.zeros(np.size(CubDrawer.xPoint, 0))
x5[:] = CubDrawer.xPoint[:,0]
y5 = np.zeros(np.size(CubDrawer.xPoint, 0))
y5[:] = CubDrawer.xPoint[:,1]

################## Regular Cubature
dim                 = 2
#rule                = 'Clenshaw-Curtis'
rule                = 'Regular'
distributionType    = 'Uniform'
xDistributionParam  = [0.0,1.0]
nPointPerDim        = 7
xRule               = [rule] * dim
xDistributionType   = [distributionType] * dim
xxDistributionParam = [xDistributionParam] * dim

#RCubDrawer = CubatureSampling(None, rule=rule, distribution=distribution, nPointPerDim=nPointPerDim, verbosity = verb+1)
RCubDrawer = CubatureSampling(None, xRule=xRule, xDistributionType=xDistributionType, xxDistributionParam=xxDistributionParam, nPointPerDim=nPointPerDim, isSparse = False, verbosity = verb+1)
RCubDrawer.Draw()

x6 = np.zeros(np.size(RCubDrawer.xPoint, 0))
x6[:] = RCubDrawer.xPoint[:,0]
y6 = np.zeros(np.size(RCubDrawer.xPoint, 0))
y6[:] = RCubDrawer.xPoint[:,1]


############## Regular Quadrature
#---Regular
rule               = 'Simpson'
maxNPoint             = 8
distributionType   = 'Uniform'
xDistributionParam = [0.0, 1.0]
x7    = []
y7    = []
for nPoint in range(2, maxNPoint+1):
  #rQuadDrawer = QuadratureSampling(None, rule, distributionType, xDistributionParam, nPoint = nPoint, verbosity=verb)
  rQuadDrawer = QuadratureSampling(None, rule, distributionType, xDistributionParam, level = nPoint, verbosity=verb)
  rQuadDrawer.Draw()
  x7Temp = rQuadDrawer.xPoint.tolist()
  y7Temp = rQuadDrawer.xPoint.tolist()
  for k in range(len(y7Temp)):
    y7Temp[k] = nPoint
  x7 = x7 + x7Temp
  y7 = y7 + y7Temp

#--------------------------------------------------------------------- Plot
## MC
#fig = plt.figure(figsize=(14,10))
#plt.scatter(x1, y1, color='green', label='MC')
#plt.legend(loc='best')
#plt.show()
#
## QMC
#fig = plt.figure(figsize=(7,5))
#plt.scatter(x2, y2, color='k',   label='QMC')
#plt.xlabel(r'$\chi_T$')
#plt.ylabel(r'$\chi_{A_{12}}$')
#plt.xlim([0.0,1.0])
#plt.ylim([-1.5,1.5])
##plt.legend(loc='best')
#plt.show()
#
## RQMC
#fig = plt.figure(figsize=(14,10))
#plt.scatter(x3, y3, label='RQMC')
#plt.legend(loc='best')
#plt.show()
#
# Quadrature
fig = plt.figure(figsize=(14,10))
plt.scatter(x4, y4, label='Quadrature')
plt.legend(loc='best')
plt.show()

## Cubature
#fig = plt.figure(figsize=(14,10))
#plt.scatter(x5, y5, label='Cubature')
#plt.legend(loc='best')
#plt.show()

## Regular Cubature
#fig = plt.figure(figsize=(14,10))
#plt.scatter(x6, y6, label='Regular Cubature')
#plt.legend(loc='best')
#plt.show()

## Regular Cubature
#fig = plt.figure(figsize=(14,10))
#plt.scatter(x7, y7, label='Regular Quadrature')
#plt.legend(loc='best')
#plt.show()
