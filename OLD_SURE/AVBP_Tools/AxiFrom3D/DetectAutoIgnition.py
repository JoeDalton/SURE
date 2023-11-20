import os
import numpy as np
import matplotlib
#matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import AVBP_Tools.AxiFrom3D.InterpolateFromAVBP as ci
import common_IO as cio
import math
from scipy.interpolate import RectBivariateSpline as i2d
from Utils.I_O               import ProgBar

def DetectAutoIgnition(field, threshold, xq, rq):
  # Thanks Corentin Grimaldi for the simplest of solutions !
  # It is possible to add an local interpolation around this estimation, but it is overkill if the interpolation mesh is already fine enough
  nx = len(xq)
  nr = len(rq)
  for i in range(nx):
    temp = field[:,i]
    for j in range(nr):
      if field[j,i] >= threshold:
        h = xq[i]
        return h
      
#def OldBuggyDetectAutoIgnition(field, threshold, xq, rq):
#  bContour = plt.contour(rq, xq, np.flipud(np.rot90(field, 1)), [threshold], colors=['w'])
#  h = min(bContour.allsegs[0][0][:, 1])
#  return h

def ComputeAverageField(lfield, ltav):
  av_field = np.zeros_like(lfield[0])
  tot_time = 0
  for idxf in range(len(lfield)):
    av_field += lfield[idxf] * ltav[idxf]
    tot_time += ltav[idxf]
  av_field /= tot_time
  return av_field
  
def BootStrapAutoIgnit(lfield, ltav, nReSample, threshold, xq, rq, smoothed=False, outputDist=True):
  data_mean = np.zeros(nReSample)
  lsize = len(lfield)
  sizeReSample = lsize
  xfield = np.array(lfield, dtype=object)
  xtav = np.array(ltav)

  xIdx = np.arange(lsize, dtype=int)
  myProgBar = ProgBar('Computing bootstrap estimate...', nReSample)
  for i in range(nReSample):  # Number of bootstrap estimates
    sample_i = np.random.choice(xIdx, size=sizeReSample, replace=True)
    field_i = ComputeAverageField(xfield[sample_i], xtav[sample_i])
    data_mean[i] = DetectAutoIgnition(field_i, threshold, xq, rq)
    myProgBar.Update()

  myProgBar.Terminate()
  if smoothed:
    standardDev = 1.0 / math.sqrt(nReSample)  # Found in wikipedia article for the smoothed bootstrap. No other source, though...
    smoothing = np.random.normal(loc=0.0, scale=standardDev, size=nReSample)
    data_mean += smoothing

  if outputDist:
    return data_mean
  else:
    return np.mean(data_mean)

