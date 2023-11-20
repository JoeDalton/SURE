#!/usr/bin/env python2.7
import argparse
import Link
import matplotlib.pyplot as plt
from AVBP_Tools.PostProc_ProdZv import *
import numpy as np
import Utils.PrettyPlots as pp


parser = argparse.ArgumentParser(description="Generates the interpolation file for an axisymetric mean of an AVBP solution")
parser.add_argument('--wd'    , nargs='+', type=str , help = "List of workdirs")
parser.add_argument('--interp', nargs=1, type=str,    help = "relative path to read the interpolation file")
parser.add_argument('--prefix', nargs=1, type=str,    help = "prefix to add to the files")
parser.add_argument('-D'      , nargs=1, type=float,  help = "Normalizing length (Nozzle diameter for ex.)")
args = parser.parse_args()

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

for i in range(len(xWorkDir)):
  PostProcProdZv(xWorkDir[i], D, interpRelPath)

