#!/usr/bin/env python2.7

import os
import sys
import Link
import numpy                                                    as np
import matplotlib.pyplot                                        as plt
import Utils.PrettyPlots                                        as pp
from Numerics.Statistics import ComputeBootstrappedEstimate     as boot








mean = 100
stdev = 10

#population = np.random.normal(loc=mean, scale=stdev, size=50000)
population = np.random.randint(0, 500, size=50000)

# We extract some samples from the population: These are the obseervations
data = population[:100]

iniMean = np.mean(data)
#bootDist= boot(data, 10000, outputDist=True)
bootDist= boot(data, 10000, smoothed=True, outputDist=True)
bootMean = np.mean(bootDist)


myKDE = pp.KDE(bootDist)
myKDE.OptimizeBandwidth(verySmooth=True)
x, y = myKDE.GetDataPlot()

print (iniMean, bootMean)

fig = plt.figure()
ax= plt.subplot(111)
ax.hist(bootDist, density=True)
ax.plot(x,y, color='k', linestyle='-')
plt.show()




#for i in range(len(xLabel)):
#  if xKDE[i]:
#    myKDE = pp.KDE(xxEval[i])
#    myKDE.OptimizeBandwidth(verySmooth=xSmKDE[i])
#    x, y = myKDE.GetDataPlot()
#    if xIsRef[i]:
#      ax.plot(x,y, label=xLabel[i], color='k', linestyle='-')
#    else:
#      ax.plot(x,y, label=xLabel[i], color='grey', linestyle=linestyles[i%len(linestyles)])
#  else:
#    ax.hist(xxEval[i], xNCol[i], density=True, alpha=0.33, label=r'' + xLabel[i])
#
#ax.set_xlabel(r'' + QoIName)
#ax.set_ylabel(r'pdf(' + QoIName + ')')
#ax.set_xlim([-8, -2])
#ax.set_title('PDF comparison')
#ax.legend(loc='best')
#pp.CorrectAxSize(ax)
##plt.show()
#plt.savefig(plotFile, bbox_inches='tight')
#
#
#prt('Distances:', 'blue', True)
#print(xDist)
#
#prt("THE END", 'green', True)
#quit()
