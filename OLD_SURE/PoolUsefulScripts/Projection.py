from Legendre import *


def BuildPCE(dim, level, xValue):


  if dim == 1:
    nCo = 3
    dummyPol = LegendreUnivariate(np.array([1]*nCo))
    xCo = np.zeros(nCo)
    for i in range(nCo):
      xCo[i] = Projection(dim, level, xValue, dummyPol, i)
    return LegendreUnivariate(xCo)
  else:
    if dim == 2:
      nCo = 6
    elif dim == 3:
      nCo = 10
    else:
      raise Exception("Max dimension is 2, min is 1")
    dummyPol = LegendreMultivariate(np.array([1]*nCo), dim)
    xCo = np.zeros(nCo)
    for i in range(nCo):
      xCo[i] = Projection(dim, level, xValue, dummyPol, i)
    return LegendreMultivariate(xCo, dim)


def Projection(dim, level, xValue, polynomial, nMonom=-1):
  xPoint, xWeight = GetPointsAndWeights(dim, level)
  result = 0
  if nMonom == -1:
    for i in range(len(xPoint)):
      result += xWeight[i] * xValue[i] * polynomial(xPoint[i])
    return result
  else:
    for i in range(len(xPoint)):
      result += xWeight[i] * xValue[i] * polynomial(xPoint[i], nMonom=nMonom)
    return result


def GetPointsAndWeights(dim, level):
  xPoint1D, xWeight1D = GetPointsAndWeights1D(level)
  if dim == 1:
    return xPoint1D, xWeight1D
  elif dim == 2:
    n = 2 * level + 1
    xPoint = np.zeros((n,n,2))
    xWeight = np.zeros((n,n))
    for i in range(n):
      for j in range(n):
        xPoint[i,j,0] = xPoint1D[i]
        xPoint[i,j,1] = xPoint1D[j]
        xWeight[i,j]  = xWeight1D[i]*xWeight1D[j]
    xWeight = xWeight.reshape(n*n)
    xPoint = xPoint.reshape((n*n,2))
  elif dim == 3:
    xPoint = np.zeros((n,n,n,3))
    xWeight = np.zeros((n,n,n))
    for i in range(n):
      for j in range(n):
        for k in range(n):
          xPoint[i,j,k,0] = xPoint1D[i]
          xPoint[i,j,k,1] = xPoint1D[j]
          xPoint[i,j,k,2] = xPoint1D[k]
          xWeight[i,j,k]  = xWeight1D[i]*xWeight1D[j]*xWeight1D[k]

    xWeight = xWeight.reshape(n*n*n)
    xPoint = xPoint.reshape((n*n*n,3))
  return xPoint, xWeight


def GetPointsAndWeights1D(level):
  nPoint = 2 * level + 1
  xMin = -1.0
  xMax = 1.0
  step = (xMax-xMin) / (nPoint-1)
  xWeight = np.zeros(nPoint)
  xPoint = np.zeros(nPoint)
  for i in range(nPoint):
    # This should work
    xi = xMin + i * step
    if i == 0 or i == nPoint -1:
      wi = (xMax - xMin) / (level * 1.0)
    elif i%2 == 1: 
      wi = 4.0 * (xMax - xMin) / (level * 1.0)
    else:
      wi = 2.0 * (xMax - xMin) / (level * 1.0)
    xWeight[i] = wi / 12.0
    xPoint[i] = xi
  return xPoint, xWeight


