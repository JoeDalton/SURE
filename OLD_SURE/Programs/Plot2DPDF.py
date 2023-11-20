#!/usr/bin/env python2.7

import os
import sys
import h5py
import Link
import argparse
import numpy                                                    as np
import openturns                                                as ot
import matplotlib.pyplot                                        as plt
import mpl_toolkits.mplot3d
import Utils.PrettyPlots                                        as pp
from Utils.I_O                      import prt
from Utils.I_O                      import DumpObject
from Utils.I_O                      import InputHandler         as IH
from Utils.I_O                      import GetObjectType
from scipy.stats                    import wasserstein_distance as wdist
from scipy.stats                    import gaussian_kde as kde
from Sampling.MC                    import MC
from Numerics.Algebra               import norm2, norminf
from Experiment.Experiment          import Experiment
from SurrogateModelling.PCE         import *
from SurrogateModelling.OneDSpline  import *
from SurrogateModelling.Kriging  import *


# ----- Parsing arguments -----
parser      = argparse.ArgumentParser(description="PDFComparator")
parser.add_argument('-f',	          required=False,	type=str,   help = "Relative path of data file")
args        = parser.parse_args()

dataPath  = 'Data/QoI'
inputPath = 'Data/Inputs'
myFile = args.f

datFile = h5py.File(myFile, 'r')
xEval = datFile[dataPath][:]
xT    = datFile[inputPath][:,31]
datFile.close()
temp = []
for dat in xEval:
  if dat != 0.0:
    temp.append(dat)
xEval = np.log10(np.array(temp))

#xmin = xT.min()
#xmax = xT.max()
#ymin = xEval.min()
#ymax = xEval.max()
#xx, yy = np.mgrid[xmin:xmax:100, ymin:ymax:100]
#
#positions = np.vstack([xx.ravel(), yy.ravel()])
#values = np.vstack([xT,xEval])
#kernel=kde(values)
#f = np.reshape(kernel(positions).T, xx.shape)


#fig = plt.figure()
#ax = plt.axes(projection='3d')
#surf = ax.plot_surface(xx, yy, f, cmap='coolwarm', edgecolor='none')
#ax.set_xlabel('T')
#ax.set_ylabel('log10(tau)')
#ax.set_zlabel('pdf')
#fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.view_init(60,35)
#plt.show()

fig = plt.figure()
ax = fig.gca()
ax.scatter(xT, xEval)
plt.show()




quit()
