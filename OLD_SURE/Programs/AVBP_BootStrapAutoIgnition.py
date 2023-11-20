#!/usr/bin/env python2.7
import argparse
import Link
import matplotlib.pyplot as plt
from AVBP_Tools.AxiFrom3D.DetectAutoIgnition import BootStrapAutoIgnit
import numpy as np
import Utils.PrettyPlots as pp
import AVBP_Tools.AxiFrom3D.common_IO as cio
import AVBP_Tools.AxiFrom3D.InterpolateFromAVBP as ci
from Utils.I_O                import prt
import h5py

parser = argparse.ArgumentParser(description="Bootstraps the lift-off height of an avbo flame")
parser.add_argument('-f', '--files', required=True, nargs='+', type=str , help = "List of decomposed means")
parser.add_argument('--var'   ,type=str , help = "Threshold variable, default is c_out")
parser.add_argument('--interp', nargs=1, type=str,    help = "relative path to read the interpolation file")
parser.add_argument('-t', '--threshold', type=float,  help = "Threshold, default is 0.1")
parser.add_argument('-D'      , type=float,  help = "Normalizing length (Nozzle diameter for ex.)")
parser.add_argument('-N'      , type=int,  help = "Number of bootstrap samples. Default is 1000")
parser.add_argument('--smoothed', type=bool,  help = "Smoothing of bootstrap estimate. Default is false")
args = parser.parse_args()

lFile = args.files
if args.interp is not None:
  interpRelPath = args.interp
else:
  interpRelPath='./interp.npz'
if args.threshold is not None:
  thresh = args.threshold
else:
  thresh = 0.0006
if args.var is not None:
  var = args.var
else:
  var = 'OH'
if args.D is not None:
  D = args.D
else:
  D = 0.00457
if args.N is not None:
  N = args.N
else:
  N = 1000
if args.smoothed is not None:
  sm = args.smoothed
else:
  sm = False

prt('Reading interpolation file...', 'green', 1)
interp_base = cio.load_interp_base(interpRelPath)
wts = interp_base['wts']
vtx = interp_base['vtx']
# Read bounding box of interpolation mesh from header
xlim = np.array(interp_base['header'][1][3])
rlim = np.array(interp_base['header'][1][4])
# Get mesh resolution from header
resx = interp_base['header'][1][0]
resr = interp_base['header'][1][1]
rest = interp_base['header'][1][2]
# Generate interpolation mesh-vectors
xq = np.linspace(xlim[0], xlim[1], resx)
rq = np.linspace(rlim[0], rlim[1], resr)





lfield  = []
ltav    = []
prt('Reading solution files...', 'green', 1)
for f in lFile:
  print(f)
  hr = cio.read_avbp_solution(f, 'Average/' + var)
  print "Interpolating"
  vq = ci.interpolate(hr, vtx, wts)

  print "axi mean"
  xQtemp = []
  for it in range(rest):
    imin = it * resx * resr
    imax = (it + 1) * resx * resr
    xQtemp.append(vq[imin:imax])
  xQtemp = np.array(xQtemp)
  mean = np.zeros(resx * resr)

  for ip in range(resx * resr):
    mean[ip] = np.sum(xQtemp[:, ip]) / rest
  Q = np.reshape(mean, (resr, resx))
  Q[np.isnan(Q)] = 0
  lfield.append(Q)
  fh5 = h5py.File(f,'r')
  tav= float(fh5.get('Parameters/t_av')[0])
  ltav.append(tav)


prt('Computing bootstrap estimate...', 'green', 1)
xMean = BootStrapAutoIgnit(lfield, ltav, N, thresh, xq, rq, smoothed=sm)
print("Mean = " + str(np.mean(xMean)/D))
print("Std = " + str(np.std(xMean)/D))
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)
ax.hist(xMean/D)
plt.show()
