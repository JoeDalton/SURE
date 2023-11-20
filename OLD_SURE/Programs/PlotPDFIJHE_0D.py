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
parser.add_argument('-f',nargs='+',	          required=False,	type=str,   help = "Relative path of data file")
args        = parser.parse_args()

myPath = 'Data/QoI'
myFile = args.f




if len(myFile) == 1:
  datFile = h5py.File(myFile[0], 'r')
  xEval = datFile[myPath][:]
  datFile.close()
  temp = []
  for dat in xEval:
    if dat > 0.0:
      temp.append(dat)
  xEval = np.log10(np.array(temp))
  
  fig = plt.figure(figsize=(7,5))
  ax = fig.add_subplot(111)
  
  
  
  prt('Plotting...', 'green', True)
  
  myKDE = pp.KDE(xEval)
  myKDE.OptimizeBandwidth(verySmooth=True)
  x, y = myKDE.GetDataPlot()
  ax.plot(x,y, color='k', linestyle='-')
  
  ax.set_xlabel(r'log$_{10}$($\tau$)')
  ax.set_ylabel(r'pdf(log$_{10}$($\tau$))')
  #ax.set_xlim([-3.5, -1])
  #ax.set_title('PDF comparison')
  #ax.legend(loc='best')
  pp.CorrectAxSize(ax)
  #plt.show()
  plotFile = 'pdf0D.pdf'
  plt.savefig(plotFile, bbox_inches='tight')

else:
  samplenum = np.zeros(len(myFile), dtype=int)
  for i in range(len(myFile)):
    samplenum[i] = int(myFile[i].split("_")[1].split("/")[0])
  idxList = np.argsort(samplenum)
  myFile = np.array(myFile) 
  myFile = myFile[idxList]
  

  #myFile.sort()

  xxEval = []
  for i in range(len(myFile)):
    datFile = h5py.File(myFile[i], 'r')
    xEval = datFile[myPath][:]
    datFile.close()
    temp = []
    for dat in xEval:
      if dat != 0.0:
        temp.append(dat)
    xEval = np.log10(np.array(temp))
    xxEval.append(xEval)
  xxEval = np.array(xxEval).T
  
  prt('Plotting...', 'green', True)
  fig = plt.figure(figsize=(30,6))
  ax = fig.add_subplot(111)
  ax.set_xlabel(r'$T_{ini}$ [K]')
  ax.set_ylabel(r'log$_{10}$($\tau$) [-]')
  #DEBUG
  #precursorList=np.linspace(990, 1050, 40)
  precursorList=np.linspace(850, 1000, 15)
  
  v1 = ax.violinplot(xxEval, points=1000,positions=precursorList,showmeans=False, showextrema=False, showmedians=False, widths=10)
  for b in v1['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further left than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
    b.set_facecolor('#D43F3A')
    b.set_edgecolor('black')
    b.set_alpha(1)


  plt.title("With Konnov 2008")
  plt.show()
  #plt.savefig(plotFile, bbox_inches='tight')


quit()
