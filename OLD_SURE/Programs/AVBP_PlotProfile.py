#!/usr/bin/env python2.7
import argparse
import Link
import matplotlib.pyplot as plt
from AVBP_Tools.AxiFrom3D.PlotProfile import PlotProfile
import numpy as np
import Utils.PrettyPlots as pp
from Utils.I_O import ErrorMsg, LoadCSV 

parser = argparse.ArgumentParser(description="Generates the interpolation file for an axisymetric mean of an AVBP solution")
parser.add_argument('--wd'    , nargs='+', type=str , help = "List of workdirs")
parser.add_argument('--var'   , nargs='+', required=True, type=str , help = "List of variables")
parser.add_argument('--label' , nargs='+', type=str , help = "List of labels corresponding to the files in each workdir")
parser.add_argument('--xNX'   , nargs='+', type=float,help = "List of axial positions for radial cuts. If not provided, an axial plot is done")
parser.add_argument('--interp', type=str,    help = "relative path to read the interpolation file")
parser.add_argument('--prefix', nargs=1, type=str,    help = "prefix to add to the files")
parser.add_argument('-D'      , nargs=1, type=float,  help = "Normalizing length (Nozzle diameter for ex.)")
parser.add_argument('--csv'   , nargs='+', type=str,    help = "profile in csv file to plot")
parser.add_argument('--csvlabel', nargs='+', type=str,  help = "label of profiles in csv")
args = parser.parse_args()

xNX = args.xNX
xVar = args.var
if args.xNX is not None:
  xNX    = args.xNX
  isAxial = False
else:
  xNX = [1]
  isAxial = True
if args.wd is not None:
  xWorkDir    = args.wd
else:
  xWorkDir = ['.']
if args.csv is not None:
  xcsv   = args.csv
else:
  xcsv    = []
if args.label is not None:
  xLabel   = args.label
else:
  xLabel    = xWorkDir
if args.csvlabel is not None:
  xcsvLabel   = args.csvlabel
else:
  xcsvLabel    = xcsv
if args.interp is not None:
  interpRelPath = args.interp
else:
  interpRelPath='interp.npz'
if args.prefix is not None:
  prefix = str(args.prefix[0]) + '_'
else:
  prefix=''
if args.D is not None:
  D = args.D
else:
  D = 0.00457

assert len(xWorkDir)==len(xLabel), "The lists of workdirs and labels must be of same length."


rstep = 2.0 # space between ticks on the abscissa of the plot

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
    print('Please check that that the y label corresponds to what you expect')
    yLab = r'$' + var + '$ [-]'
  fig = plt.figure(figsize=(8.0,8.0*len(xNX)))
  nPlot = len(xNX)
  xAx = []
  for i in range (nPlot):
    xAx.append(plt.subplot(nPlot, 1, i+1))

  for i in range(len(xWorkDir)):
    PlotProfile(xWorkDir[i], D, var, xLabel[i], xAx, xNX, interpRelPath, axial=isAxial)

  print "Prepare plot and save figure"
  #rc = {"figure.subplot.left"    : 0.1,
  #             "figure.subplot.right"   : 0.9,
  #             "figure.subplot.bottom"  : 0.1,
  #             "figure.subplot.top"     : 0.9 }
  if isAxial:
    maxX = 50.0*D
    xticks = np.arange(0.0, maxX, rstep * D * 2.5)
  else:
    maxX = 10.0*D
    xticks = np.arange(0.0, maxX, rstep * D)
  xtickslabels = np.round(xticks/D)
  for i in range(len(xAx)):
    xAx[i].set_aspect('auto')
    xAx[i].set_xlim([0.0, maxX])
    xAx[i].set_xticks(xticks)
    xAx[i].set_xticklabels(xtickslabels)
    if isAxial:
      xAx[i].set_xlabel('x/D')
    else:
      xAx[i].set_xlabel('r/D')
      xAx[i].set_title('X = ' + str(xNX[i]) + ' D')
    xAx[i].set_ylabel(yLab)
    pp.CorrectAxSize(xAx[i])
  pp.AdjustSubplotSpacings(doubleVert=True)

  if len(xNX) == 1 and len(xVar) == 1:
    for icsv in range(len(xcsv)):
      x,y = LoadCSV(xcsv[icsv])
      if icsv == 0:
        xAx[0].scatter(x*D,y, label=xcsvLabel[icsv], edgecolors='k',facecolors='none', marker='o')
        print(x[-1]*D)
      else:
        xAx[0].plot(x,y, label=xcsvLabel[icsv])
        print(x[-1])
    xAx[0].legend(fontsize=12)
  else:
    ErrorMsg("I can't add csv profiles if several variables or axial locations are selected", "AVBP_PlotProfile.py")


  plt.savefig('Profiles_' + prefix + var + '.pdf', bbox_inches="tight", pad_inches = 0)
  fig = None
