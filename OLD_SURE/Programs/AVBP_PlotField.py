#!/usr/bin/env python2.7
import argparse
import Link
import matplotlib.pyplot as plt
from AVBP_Tools.AxiFrom3D.PlotField import PlotField
import numpy as np
import Utils.PrettyPlots as pp


parser = argparse.ArgumentParser(description="Generates the interpolation file for an axisymetric mean of an AVBP solution")
parser.add_argument('--wd'    , nargs='+', type=str , help = "List of workdirs")
parser.add_argument('--var'   , nargs=1, required=True, type=str , help = "List of variables")
parser.add_argument('--interp', type=str,    help = "relative path to read the interpolation file")
parser.add_argument('--prefix', nargs=1, type=str,    help = "prefix to add to the files")
parser.add_argument('-D'      , nargs=1, type=float,  help = "Normalizing length (Nozzle diameter for ex.)")
args = parser.parse_args()

xVar = args.var
if args.wd is not None:
  xWorkDir    = args.wd
else:
  xWorkDir = ['.']
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

var = xVar[0]
for i in range(len(xWorkDir)):
  print(xWorkDir[i])
  if var == "c_out":
    PlotField(xWorkDir[i], D, var, interpRelPath, xIso=[0.1], bottomThresh=0.1)
  elif var == "OH":
    PlotField(xWorkDir[i], D, var, interpRelPath, xIso=[], bottomThresh=0.0006)
    #PlotField(xWorkDir[i], D, var, interpRelPath, xIso=[], bottomThresh=None)
  else:
    PlotField(xWorkDir[i], D, var, interpRelPath, xIso=[], bottomThresh=None)

