#!/usr/bin/env python2.7
import argparse
import Link
import matplotlib.pyplot as plt
from AVBP_Tools.AxiFrom3D.GetField import GetField
import numpy as np
import Utils.PrettyPlots as pp
import h5py

parser = argparse.ArgumentParser(description="Stores an AVBP solution into a SURE result file")
parser.add_argument('--res'   , required=True, type=str , help = "Result file")
parser.add_argument('--wd'    , nargs='+', type=str , help = "List of workdirs to look for solutions")
parser.add_argument('--var'   , required=True, type=str , help = "Variable to look for")
parser.add_argument('--interp', type=str,    help = "relative path to read the interpolation file")
args = parser.parse_args()

res = args.res
var = args.var
if args.wd is not None:
  xWorkDir    = args.wd
else:
  xWorkDir = ['.']
if args.interp is not None:
  interpRelPath = args.interp
else:
  interpRelPath='interp.npz'


resFile  = h5py.File(res, 'r+')
for i in range(len(xWorkDir)):
    mean, resx, resr, xlim, rlim = GetField(xWorkDir[i], var, interpRelPath)
    print "Exporting...."
    resFile.create_dataset('Fields/' + xWorkDir[i] + '/' + var + "/Array", data = mean)
    resFile.create_dataset('Fields/' + xWorkDir[i] + '/' + var + "/resx", data = resx)
    resFile.create_dataset('Fields/' + xWorkDir[i] + '/' + var + "/resr", data = resr)
    resFile.create_dataset('Fields/' + xWorkDir[i] + '/' + var + "/xlim", data = xlim)
    resFile.create_dataset('Fields/' + xWorkDir[i] + '/' + var + "/rlim", data = rlim)

