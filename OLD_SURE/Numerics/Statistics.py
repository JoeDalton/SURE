from Utils.I_O import prt
import math
import sys
import numpy as np
from scipy.stats import energy_distance
from Numerics.Algebra import norm2


def ComputeDistanceFromSamples(xSample1, xSample2):
  return energy_distance(xSample1, xSample2)


def ComputeNormalizedDistanceFromSamples(xRef, xTested):
  norm      = math.sqrt(2) * norm2(xRef)
  distance  = ComputeDistanceFromSamples(xRef, xTested)
  return distance/norm

def ComputeBootstrappedEstimate(data, nReSample, sizeReSample=None, smoothed=False, outputDist=False):
  data_mean = np.zeros(nReSample)
  if sizeReSample is None:
    sizeReSample = len(data)

  for i in range(nReSample): #Number of bootstrap estimates
    sample_n = np.random.choice(data, size=sizeReSample, replace=True)
    data_mean[i] = sample_n.mean()

  if smoothed:
    standardDev = 1.0 / math.sqrt(nReSample) # Found in wikipedia article for the smoothed bootstrap. No other source, though...
    smoothing = np.random.normal(loc=0.0, scale=standardDev, size=nReSample)
    data_mean += smoothing

  if outputDist:
    return data_mean
  else:
    return np.mean(data_mean)

def WeightedMean(data, weights):
  assert (len(data)==len(weights)), "Incompatible lengths for data and weights"
  sumwgt = np.sum(weights)
  assert(sumwgt!=0.0), "Sum of weights cannot be 0"
  result = 0.0
  for i in range(len(data)):
    result += data[i] * weights[i]
  result/=sumwgt
  return result

def WeightedStdDev(data, weights):
  assert (len(data)==len(weights)), "Incompatible lengths for data and weights"
  sumwgt  = np.sum(weights)
  assert(sumwgt!=0.0), "Sum of weights cannot be 0" # Verifies indirectly that m != 0
  m       = np.count_nonzero(weights)
  mean    = WeightedMean(data, weights)
  result  = 0.0
  for i in range(len(data)):
    result += (data[i]-mean) * (data[i]-mean) * weights[i]
  denominator = (m-1) * sumwgt / m
  result/=denominator
  result = math.sqrt(result)
  return result
