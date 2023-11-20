from Utils.I_O import prt, ErrorMsg, WarningMsg
import math
import sys
import numpy as np
from numpy.linalg import matrix_rank as rank
import copy

fileName = 'Numerics/Algebra.py'

def CheckVectorsEquality(a, b, rtol=0.0):
  # a and b must be lists or np.arrays of same size or scalars
  # TODO: Do this cleanly : "try, except" may let string types pass, for example, and it causes problems with sys.exit() which is an exception
  areEqual = True
  try:
    if len(a) != len(b):
      ErrorMsg('Vectors a and b must have the same size', fileName)
    for i in range(len(a)):
      if abs(a[i]-b[i]) > rtol:
        areEqual = False
        break
  except:
    if abs(a-b) > rtol:
      areEqual = False
  return areEqual


def ComputeActiveVectors(xxActiveVariables, xBasis):
  #Iterative routine to compute the active vectors from a list of active variables at each stage and their corresponding bases
  xxAV    = copy.deepcopy(xxActiveVariables)
  xB      = copy.deepcopy(xBasis)
  nSample = len(xxAV[0])
  xVect   = []
  iteNum  = 0
  maxIte  = len(xB)

  #inside the loop:
  while len(xB)>0:
    dim       = xB[-1].shape[0]
    initVect  = np.array([1] + [0]*(dim-1))
    vect      = np.matmul(xB[-1].T, initVect)
    temp1     = [vect]*nSample
    for idxS in range(nSample):
      temp1[idxS]   = temp1[idxS]*xxAV[-1][idxS]
      if iteNum !=0:
        temp2 = np.array([0]+list(xVect[idxS]))
        temp1[idxS] += temp2
    xVect = copy.deepcopy(temp1)

    #Prepare following ite
    iteNum += 1
    del(xB[-1])
    del(xxAV[-1])

  return np.array(xVect)


def CompleteBasis(partialBasis):
  # Partial basis must be of shape np.array([b1, ..., bc]) avec bi = [bi1, ..., bin], n > c
  # Id est np.shape(partialBasis) = (c,n)
  # (We already got some basis vectors, but we need n-c new vectors to complete the basis)
  c,n = np.shape(partialBasis)
  # A priori remaining vectors:
  xNewVec = []
  for i in range(n-c):
    newVec = [0]*i + [1] + [0]*(n-c-i)
    xNewVec.append(newVec)
  # Defining the new, completed basis:
  completeBasis = np.zeros((n,n))
  completeBasis[:c,:] = partialBasis[:,:]
  for i in range(n-c):
    completeBasis[c+i,:] = xNewVec[i][:]
  basisRank = rank(completeBasis)
  if basisRank != n:
    ErrorMsg("Dear future self, this vector family is ill-conditionned. You must implement a random draw of new basis vectors.", fileName)
  return completeBasis


def Orthonormalize(xVect, diagnostics=False, threshold=1e-08, maxIte=3):
  # xVect must be a list of non-zero numpy vectors of same length
  score       = 1.0
  counter     = 0
  result      = []
  temp1       = list(xVect) # This does a copy of xVect
  temp2       = []
  while counter < maxIte:
    temp2, score = _OrthonormalizeOnce(temp1)
    if score < threshold:
      break
    temp1 = list(temp2) # This does a copy of temp2 in temp1
    counter += 1
  
  result = list(temp2)
  if counter >= maxIte:
    WarningMsg("Max iteration of orthonormalization reached", fileName,1)

  if diagnostics:
    if len(result) > 1:
      prt('', 'red', True)
      prt('------------------------- Diagnostics', 'red', True)
      prt('Norms of the basis vectors:', 'None', True)
      for vIdx in range(len(result)):
        prt(norm2(result[vIdx]), 'None', True)
      print('Cos of two consecutive vectors:')
      for vIdx in range(len(result)-1):
        prt(result[vIdx].dot(result[vIdx+1]), 'None', True)
      prt('-------------------------', 'red', True)
  
  return np.array(result)

def _OrthonormalizeOnce(xVect):
  # xVect must be a list of non-zero numpy vectors of same length
  result = []
  score = 0.0
  if len(xVect)==0:
    return result, score
  else:
    result.append(xVect[0] / norm2(xVect[0]))
    for i in range(1, len(xVect)): # Does nothing if ther is only one vector in the list
      temp    = np.zeros_like(xVect[i])
      temp[:] = xVect[i][:]
      for vect in result:
        temp -= VectorProjection(temp, vect)
      temp = temp/norm2(temp)
      result.append(temp)

    if len(result) > 1:
      #TODO : Check instead the cos of every pair of vectors of the basis
      scoreList = []
      for vIdx in range(len(result)-1):
        scoreList.append(result[vIdx].dot(result[vIdx+1]))
      score = norminf(scoreList)

    return result, score 


def CheckOrthonormalMatrix(matrix, tol=1e-8):
  temp = np.matmul(matrix, matrix.T)
  if abs(temp-1.0)<tol:
    return True
  else:
    return False


def VectorProjection(projected, target):
  # projected and target must be numpy vectors of equal length
  #TODO: Add a check of equal length of the vectors
  result = np.zeros_like(target)
  result[:] = (projected.dot(target) / target.dot(target)) * target[:]
  return result


def ChangeBasis(vector, transferMatrix):
  return np.matmul(transferMatrix, vector)



def norm1(vect):
  return np.sum(np.absolute(np.array(vect)))

def norminf(vect):
  return np.amax(np.absolute(np.array(vect)))

def norm2(vect):
  return math.sqrt(np.array(vect).dot(np.array(vect)))


def NRMSE(refVect, vect):
  # Normalized root mean squared error
  return norm2(refVect-vect) / norm2(refVect)

def NMAE(refVect, vect):
  # Normalized maximum absolute error
  return norminf(refVect-vect) / norminf(refVect)
