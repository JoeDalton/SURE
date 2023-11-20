#!/usr/bin/env python2.7
import argparse
import Link
import matplotlib.pyplot as plt
from AVBP_Tools.AxiFrom3D.PlotProfile import MasterPlotProfile
import numpy as np
import Utils.PrettyPlots as pp
from Utils.I_O import ErrorMsg, LoadCSV 

parser = argparse.ArgumentParser(description="Generates the interpolation file for an axisymetric mean of an AVBP solution")
parser.add_argument('--wd'    , nargs='+', type=str , help = "List of workdirs")
parser.add_argument('--refd'  , type=str , help = "Directory where the reference is")
parser.add_argument('--label' , nargs='+', type=str , help = "List of labels corresponding to the files in each workdir")
parser.add_argument('--csvdir', type=str)
#parser.add_argument('--csv'   , nargs='+', type=str,    help = "profile in csv file to plot")
#parser.add_argument('--csvlabel', nargs='+', type=str,  help = "label of profiles in csv")
args = parser.parse_args()

if args.wd is not None:
  xWorkDir    = args.wd
else:
  xWorkDir = ['.']

if args.refd is not None:
  refDir    = args.refd
else:
  refDir = '../REF'

if args.csvdir is not None:
  csvDir    = args.csvdir
else:
  csvDir = '../../ExpProfiles/Cabra2004Stripped'
  #csvDir = '../ExpProfiles/Cabra2004Stripped'
  #csvDir = '../../../ExpProfiles/Cabra2004Stripped'

#if args.csv is not None:
#  xcsv   = args.csv
#else:
#  xcsv    = []
if args.label is not None:
  xLabel   = args.label
else:
  xLabel    = xWorkDir

#xWorkDir = [refDir] + xWorkDir
#xLabel = ['Reference'] + xLabel

xWorkDir = xWorkDir
xLabel = xLabel
#if args.csvlabel is not None:
#  xcsvLabel   = args.csvlabel
#else:
#  xcsvLabel    = xcsv




D = 0.00457

assert len(xWorkDir)==len(xLabel), "The lists of workdirs and labels must be of same length."


rstep = 2.0 # space between ticks on the abscissa of the plot



xVar  = ['T', 'H2', 'O2', 'H2O', 'OH', 'z_out'] #'u']

xNX   = [1,9,11,14,26]


for var in xVar:
  if var=='T':
    yLab = r'Temperature [K]'
  elif var == 'Yzv':
    yLab = r'Zv [-]'
  elif var == 'z_out':
    yLab = r'Z [-]'
  elif var=='c_out':
    yLab = r'C [-]'
  else:
    print('Please check that the y label corresponds to what you expect')
    yLab = r'$' + var + '$ [-]'
  fig = plt.figure(figsize=(8.0, 7.0 * ( len(xNX) + 1 )) )
  nPlot = len(xNX) + 1
  xAx = []
  for i in range (nPlot):
    xAx.append(plt.subplot(nPlot, 1, i+1))

  for i in range(len(xWorkDir)):
    MasterPlotProfile(xWorkDir[i], D, var, xLabel[i], xAx, xNX, "interp.npz")

  print "Prepare plot and save figure"
  #rc = {"figure.subplot.left"    : 0.1,
  #             "figure.subplot.right"   : 0.9,
  #             "figure.subplot.bottom"  : 0.1,
  #             "figure.subplot.top"     : 0.9 }
  for i in range(len(xAx)):
    if i == 0:
      maxX = 50.0*D
      xticks = np.arange(0.0, maxX, rstep * D * 2.5)
      xAx[i].set_xlabel('x/D')
      xAx[i].set_title('Axial')
    else:
      maxX = 10.0*D
      xticks = np.arange(0.0, maxX, rstep * D)
      xAx[i].set_xlabel('r/D')
      xAx[i].set_title('X = ' + str(xNX[i-1]) + ' D')
    xtickslabels = np.round(xticks/D)
    xAx[i].set_aspect('auto')
    xAx[i].set_xlim([0.0, maxX])
    xAx[i].set_xticks(xticks)
    xAx[i].set_xticklabels(xtickslabels)
    xAx[i].set_ylabel(yLab)
    pp.CorrectAxSize(xAx[i])
  pp.AdjustSubplotSpacings(doubleVert=True)

  # Add exp data from csvs
  if var=="z_out":
    csvPrefx = csvDir + '/' + 'Z' + '_'
  else:
    csvPrefx = csvDir + '/' + var + '_'
  csvSufx = '.csv'
  
  # Axial
  csvName = csvPrefx + 'Axial' + csvSufx
  x,y = LoadCSV(csvName)
  xAx[0].scatter(x*D,y, label='Exp', edgecolors='k',facecolors='none', marker='o')
  xAx[0].legend(fontsize=12)

  # Radial
  for i in range(len(xNX)):
    csvName = csvPrefx + 'z:d=' + str(xNX[i]) + csvSufx
    x,y = LoadCSV(csvName)
    xAx[i+1].scatter(np.abs(x/1000.0),y, label='Exp', edgecolors='k',facecolors='none', marker='o')
    xAx[i+1].legend(fontsize=12)

  plt.savefig('Master_Profiles_X_' + var + '.pdf', bbox_inches="tight", pad_inches = 0)
  fig = None
