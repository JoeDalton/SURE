#!/usr/bin/env python2.7

import os
import sys
import h5py
import Link
import argparse
import numpy                                                    as np
import openturns                                                as ot
import matplotlib.pyplot                                        as plt
import Utils.PrettyPlots                                        as pp
from Utils.I_O                      import prt
from Utils.I_O                      import DumpObject
from Utils.I_O                      import InputHandler         as IH
from Utils.I_O                      import GetObjectType
from scipy.stats                    import wasserstein_distance as wdist
from Sampling.MC                    import MC
from Numerics.Algebra               import norm2, norminf
from Experiment.Experiment          import Experiment
from SurrogateModelling.PCE         import *
from SurrogateModelling.OneDSpline  import *
from SurrogateModelling.Kriging  import *

# ----- Parsing arguments -----
parser      = argparse.ArgumentParser(description="PDFComparator")
parser.add_argument('-i',	          required=False,	type=str,   help = "Relative path of input file")
args        = parser.parse_args()

# ----- Initializing objects from input file -----
if args.i is not None:
  inputFile		= args.i
  INPUT = IH(inputFile, modify=False)
else:
  inputFile = 'PDFComparator.input'
  try:
    INPUT = IH(inputFile, modify=False)
  except:
    prt('Input File does not exist', 'None', True)
    prt('  --> Copy an example in ' + inputFile, 'None', True)
    os.system("cp " + Link.sureHome + '/InputFiles/' + inputFile + ' .')
    sys.exit()

##############################
#       Initialisation       #
##############################
xxEval  = []
xLabel  = []
xKDE    = []
xSmKDE  = []
xNCol   = []
xIsRef  = []
xDist   = []

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)
INPUT.GetBlock('PDFComparator')
keys      = INPUT.GetKeysForBlock()
plotFile  = INPUT.GetCharacterFromKey('PlotFile')
QoIName   = INPUT.GetCharacterFromKey('QoIName')
if 'BlackAndWhite' in keys:
  baw     = INPUT.GetLogicalFromKey('BlackAndWhite')
else:
  baw     = False
if 'CircleLineStyles' in keys:
  cls = INPUT.GetLogicalFromKey('CircleLineStyles')
else:
  cls = False
subBlocks = INPUT.GetSubBlockNames()


##############################
#      Loading DataSets      #
##############################
prt('Loading datasets', 'green', True)
counter = 0 #DEBUG
for sb in subBlocks:
  INPUT.GetSubBlock(sb)
  keys = INPUT.GetKeysForBlock()
  myFile  = INPUT.GetCharacterFromKey('DatFile')
  myPath  = INPUT.GetCharacterFromKey('DatPath')
  myLabel = INPUT.GetCharacterFromKey('DatLabel')
  if 'PlotKDE' in keys:
    plotKDE = INPUT.GetLogicalFromKey('PlotKDE')
  else:
    plotKDE = False
  if plotKDE:
    nCol = 0
    if 'SmoothenKDE' in keys:
      smKDE = INPUT.GetLogicalFromKey('SmoothenKDE')
    else:
      smKDE = True
  else:
    smKDE = False
    if 'BinNumber' in keys:
      nCol = INPUT.GetIntegerFromKey('BinNumber')
    else:
      nCol = 50
  if 'IsReference' in keys:
    isRef = INPUT.GetLogicalFromKey('IsReference')
  else:
    isRef = False

  datFile = h5py.File(myFile, 'r')
  if isRef:
    xEval = datFile[myPath][:]
  else:
    xEval = datFile[myPath][:]
  datFile.close()

  ##EMERGENCY
  #datFile = h5py.File(myFile, 'r')
  #xEval = datFile[myPath][:]
  #datFile.close()
  #temp = []
  #if counter<1:
  #  for dat in xEval:
  #    if dat != 0.0 and dat < 1.0:
  #      temp.append(dat)
  #  xEval = np.log10(np.array(temp))
  #else:
  #  xEval = np.log10(np.exp(xEval))
  #counter+=1

  #print(len(xEval))

  ## END EMERGENCY
  xxEval.append(xEval)
  xLabel.append(myLabel)
  xKDE.append(plotKDE)
  xSmKDE.append(smKDE)
  xNCol.append(nCol)
  xIsRef.append(isRef)

##############################
#    Checking ref unicity    #
##############################
prt('Checking reference unicity', 'green', True)
nRef  = 0
idRef = -1
for i in range(len(xIsRef)):
  if xIsRef[i]:
    nRef += 1
    idRef = i
if nRef == 0:
  ErrorMsg('No reference PDF specified', 'PDFComparator.py')
if nRef > 1:
  ErrorMsg('Only one reference PDF can be specified', 'PDFComparator.py')

##############################
#     Computing distances    #
##############################
prt('Computing pdf distances', 'green', True)
for i in range(len(xLabel)):
  xDist.append(wdist(xxEval[i], xxEval[idRef]))


##############################
#     Computing distances    #
##############################
prt('Plotting...', 'green', True)

linestyles=['-', '--', '-.', ':']

for i in range(len(xLabel)):
  if xKDE[i]:
    myKDE = pp.KDE(xxEval[i])
    myKDE.OptimizeBandwidth(verySmooth=xSmKDE[i])
    x, y = myKDE.GetDataPlot()
    if xIsRef[i]:
      ax.plot(x,y, label=xLabel[i], color='k', linestyle='-')
    else:
      ax.plot(x,y, label=xLabel[i], color='grey', linestyle=linestyles[i%len(linestyles)])
  else:
    ax.hist(xxEval[i], xNCol[i], density=True, alpha=0.33, label=r'' + xLabel[i])

ax.set_xlabel(r'' + QoIName)
ax.set_ylabel(r'pdf(' + QoIName + ')')
#ax.set_xlim([-3.5, -1])
#ax.set_title('PDF comparison')
ax.legend(loc='best')
pp.CorrectAxSize(ax)
#plt.show()
plt.savefig(plotFile, bbox_inches='tight')


prt('Distances:', 'blue', True)
print(xDist)

prt("THE END", 'green', True)
quit()
