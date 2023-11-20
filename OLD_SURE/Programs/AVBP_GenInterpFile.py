#!/usr/bin/env python2.7
import argparse
import Link
from AVBP_Tools.AxiFrom3D.GenerateInterpolationFile import GenInterpFile
import os

parser = argparse.ArgumentParser(description="Generates the interpolation file for an axisymetric mean of an AVBP solution")
parser.add_argument('--wd'    ,          type=str , help = "Workdir")
parser.add_argument('--mesh'           , type=str,  help = "relative path to mesh file")
parser.add_argument('--interp'         , type=str,  help = "relative path to write the interpolation file")
parser.add_argument('--lims'  , nargs=6, type=float, help = "lims in the interpolation mesh")
parser.add_argument('--cart'  , action='store_true', help = "Activates cartesian mesh instead of axisym")
parser.add_argument('--slice' , action='store_true', help = "Activates slice instead of 3d")
args = parser.parse_args()

cart=args.cart
Slice=args.slice
if args.wd is not None:
  xWorkDir    = args.wd
else:
  xWorkDir = '.'
if args.mesh is not None:
  meshRelPath = args.mesh
else:
  meshRelPath="mesh.mesh.h5"
if args.interp is not None:
  interpRelPath = args.interp
else:
  interpRelPath='interp.npz'
if args.lims is not None:
  lim1 = [float(args.lims[0]), float(args.lims[1])]
  lim2 = [float(args.lims[2]), float(args.lims[3])]
  lim3 = [float(args.lims[4]), float(args.lims[5])]
else:
  #Defaults for the Cabra flame
  if cart:
    lim1 = [0.0, 0.12]
    lim2 = [-0.01371,0.01371]
    lim3 = [-0.01371, 0.01371]
  else:
    lim1 = [0.0, 0.12]
    lim2 = [0.0, 0.01371]
    lim3 = [0.0, 360.0]


for workDir in xWorkDir:
  GenInterpFile(workDir, meshRelPath=meshRelPath, interpRelPath=interpRelPath, lim1=lim1, lim2=lim2, lim3=lim3, cart=cart, Slice=Slice)

