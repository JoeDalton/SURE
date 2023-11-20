from Utils.I_O import ErrorMsg
from scipy.special import factorial
import sys

fileName = 'Numerics/Combinatorics.py'

def ComputeFactorial(n):
  return factorial(n, exact = True)

def ComputeBinomialCoefficient(n,k):
  # For any integers n, k >= 0
  #if not n.is_integer():
  #  print("ComputeBinomialCoefficient : n must be an integer")
  #  sys.exit()
  if n<0:
    ErrorMsg("ComputeBinomialCoefficient : n must be positive", fileName)
  #if not k.is_integer():
  #  print("ComputeBinomialCoefficient : k must be an integer")
  #  sys.exit()
  if k<0:
    ErrorMsg("ComputeBinomialCoefficient : k must be positive", fileName)
  if n-k<0:
    return 0.0
  x = ComputeFactorial(n) / (ComputeFactorial(k) * ComputeFactorial(n-k))
  return x
