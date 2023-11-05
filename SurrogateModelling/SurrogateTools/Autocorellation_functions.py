import numpy as np
from math import sqrt, exp

def _Prepare(X1, X2, Theta):
  # Details: Ensures the correct shape for matrix X2 and computes sizes
  #
  # inputs:
  # X1  - First input vector 
  # X2  - Second input vector
  # Theta - hyperparameter vector
  #
  # outputs:
  # x1  - First input vector
  # x2  - Second input vector, possibly reformated
  # n   - Dimension of input in x1
  # m1  - Number of vectors in X1, the initial sample set
  # m2  - Number of vectors in X2
  # theta - hyperparameter vector

  if not hasattr(Theta, '__len__'): # 1D case, when the hyperparameter vector has not been conditionned as a vector
    theta = [Theta]
  else:
    theta = Theta

  try:
    n     = np.shape(X1)[1]
  except IndexError:
    n = 1 # Case where X1 is a vector of 1D inputs
    X1 = np.array([X1]).T # We need to take the transpose of X1 (thanks numpy)
    X2 = np.array([X2]).T # Same for X2 (thanks numpy)
  m1    = np.shape(X1)[0]
  try:
    m2  = np.shape(X2)[0]
  except:
    m2  = 1
  # Case X2 = vector of 1 input with multiple parameters
  if len(np.shape(X2))==1:
    temp = np.zeros((1,m2))
    temp[0,:] = X2[:]
    X2 = temp
    m2 = 1
  # Reformat X2
  if np.size(X2) == 1:
    #x2 = np.array([x2])
    X2 = np.array([[X2]])
  return X1, X2, n, m1, m2, theta


def Matern52_matrix(X1, X2, theta):
  # Details: Create autocorrelation matrix with Matern52 kernel
  #
  # inputs:
  # x1 - First input vector 
  # x2 - Second input vector
  # theta - Hyperparameter vector
  #
  # outputs:
  # R - Autocorrelation matrix

  
  x1, x2, n, m1, m2, theta = _Prepare(X1, X2, theta)

  R = np.ones((m1,m2))
  for i in range(m1):
    for j in range(m2):
      R_val = 1.0
      for p in range(n):
        r     =   x1[i,p]-x2[j,p]
        kval  =   (sqrt(5) * abs(r)) / theta[p]
        val   =   ( 1 + kval + kval*kval/3 ) * exp(-kval)
        R_val *=  val
      R[i,j] = R_val
  return R


def Matern32_matrix(X1, X2, theta):
  # Details: Create autocorrelation matrix with Matern32 kernel
  #
  # inputs:
  # x1 - First input vector 
  # x2 - Second input vector
  # theta - Hyperparameter vector
  #
  # outputs:
  # R - Autocorrelation matrix
  
  x1, x2, n, m1, m2, theta = _Prepare(X1, X2, theta)

  R = np.ones((m1,m2))
  for i in range(m1):
    for j in range(m2):
      R_val = 1.0
      for p in range(n):
        r     =   x1[i,p]-x2[j,p]
        kval  =   (sqrt(3) * abs(r)) / theta[p]
        val   =   (1 + kval) * exp(-kval)
        R_val *=  val
      R[i,j] = R_val
  return R


def Ex_matrix(X1, X2, theta):
  # Details: Create autocorrelation matrix with exponential kernel
  #
  # inputs:
  # x1 - First input vector 
  # x2 - Second input vector
  # theta - Hyperparameter vector
  #
  # outputs:
  # R - Autocorrelation matrix
  
  x1, x2, n, m1, m2, theta = _Prepare(X1, X2, theta)

  R = np.ones((m1,m2))
  for i in range(m1):
    for j in range(m2):
      val = 0.0
      for p in range(n):
        r     =   x1[i,p]-x2[j,p]
        val   +=  abs(r) / sqrt(theta[p])
      R[i,j] = exp(-val)
  return R


def SqEx_matrix(X1, X2, theta):
  # Details: Create autocorrelation matrix with squared exponential kernel
  #
  # inputs:
  # x1 - First input vector 
  # x2 - Second input vector
  # theta - Hyperparameter vector
  #
  # outputs:
  # R - Autocorrelation matrix
  
  x1, x2, n, m1, m2, theta = _Prepare(X1, X2, theta)

  R = np.ones((m1,m2))
  for i in range(m1):
    for j in range(m2):
      val = 0.0
      for p in range(n):
        r     =   x1[i,p]-x2[j,p]
        temp  =   r / sqrt(theta[p])
        val   +=  temp * temp
      R[i,j] = exp(-val)
  return R


def CubicSpline_matrix(X1, X2, theta):
  # Details: Create autocorrelation matrix with cubic spline kernel
  #
  # inputs:
  # x1 - First input vector 
  # x2 - Second input vector
  # theta - Hyperparameter vector
  #
  # outputs:
  # R - Autocorrelation matrix 
  
  x1, x2, n, m1, m2, theta = _Prepare(X1, X2, theta)

  R = np.ones((m1,m2))
  for i in range(m1):
    for j in range(m2):
      R_val = 1.0
      for p in range(n):
        xi = abs(x1[i,p]-x2[j,p]) / theta[p]
        if (0.0 <= xi) and (xi <= 0.2):
            val = 1 - 15.0 * xi*xi + 30.0 * xi*xi*xi
        elif (0.2 < xi) and (xi < 1.0):
            temp = 1-xi
            val = 1.25 * temp*temp*temp
        elif xi >= 1.0:
            val = 0.0
        R_val *=  val
      R[i,j] = R_val
  return R


def Linear_matrix(X1, X2, theta):
  # Details: Create autocorrelation matrix with linear kernel
  #
  # inputs:
  # x1 - First input vector 
  # x2 - Second input vector
  # theta - Hyperparameter vector
  #
  # outputs:
  # R - Autocorrelation matrix
  
  x1, x2, n, m1, m2, theta = _Prepare(X1, X2, theta)

  R = np.ones((m1,m2))
  for i in range(m1):
    for j in range(m2):
      val = 1.0
      for p in range(n):
        r     =   abs(x1[i,p]-x2[j,p])
        val   *=  max( 0, 1 - (r / theta[p]))
      R[i,j] = val
  return R
