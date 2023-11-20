#!/usr/bin/env python2.7
import argparse
import Link
from AVBP_Tools.Scatterplot import Scatterplot

parser = argparse.ArgumentParser(description="Scatter plot of AVBP solution.\n + 1 var => Histogram\n + 2 vars => scatterplot (v2=f(v1))\n + 3 vars => scatterplot (v2=f(v1)) coloured by v3")
parser.add_argument('-s', '--sol'   , type=str , help = "Solution file. Default is 'current.h5'")
parser.add_argument('-m', '--mesh'  , type=str,  help = "Mesh file. Default is 'mesh.mesh.h5'")
parser.add_argument('-t', '--title' , type=str,  help = "Title of the plot. Default is 'Scatter plot'")
parser.add_argument('-v', '--var'   , nargs='+', required=True, type=str , help = "List of vars, X, (Y), (Z)")
args = parser.parse_args()

if args.sol is not None:
  sol     = args.sol
else:
  sol     = 'current.h5'

if args.mesh is not None:
  mesh    = args.mesh
else:
  mesh  = "mesh.mesh.h5"

if args.title is not None:
  title = args.title
else:
  title = "Scatter plot"

xV = args.var

if len(xV)>0 and len(xV)<4:
  Scatterplot(sol, mesh, title, xV)
else:
  print('--var has min. 1 argument and max. 3 arguments')
