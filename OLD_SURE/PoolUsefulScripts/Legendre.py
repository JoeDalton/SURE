#!/usr/bin/env python2.7

from __future__ import print_function
import time
import numpy                                      as np
import math



class Polynomial():
  """
  Univariate Legendre polynomials of total order 2
  coefs defines the monomials of the Legendre basis
  """
  
  def __init__(self): 
    raise Exception("Base class cannot be created") 
 

  def __call__(self, x, nMonom =-1):
    if nMonom == -1:
      return self.Evaluate(x)
    else:
      return self.EvaluateMonomial(x, nMonom)


  def Evaluate(self, x):
    raise Exception("Base class cannot be evaluated") 
    return None





class LegendreUnivariate(Polynomial):
  """
  Univariate Legendre polynomials of total order 2
  coefs defines the monomials of the Legendre basis
  """
  coefs = None
  Pol = None
  
  def __init__(self, Pol): 
    """
    Pol is a list of the coefficients for the monomials. It defines completely the polynomials 
    """
    self.coefs = np.zeros((3,3))
    self.coefs[0,0] = 1
    self.coefs[1,1] = 1.732051
    self.coefs[2,0] = -1.11803
    self.coefs[2,2] = 3.3541
    #self.coefs[0,0] = 1
    #i = 0
    #self.coefs[1,1] = math.sqrt( (2*i+1)*(2*i+3) )/( i+1 )
    #i = 1
    #self.coefs[2,0] = - i * math.sqrt(2*i+3) / ( (i+1) * math.sqrt(2*i-1) )
    #self.coefs[2,2] = math.sqrt( (2*i+1)*(2*i+3) )/( i+1 ) * self.coefs[1,1]
    if len(Pol) == 3:
      self.Pol = Pol
    else:
      raise Exception("3 coeffs are needed to define a 1D polynomial of order 2") 
 

  def Evaluate(self, x):
    result = 0
    for i in range(3):
      result += self.EvaluateMonomial(x, i)
      #print('==========')
      #print(x)
      #print(result)
      #print('==========')
    return result


  def EvaluateMonomial(self, x, i):
    """
    x : evaluation point
    n : index of the monomial in the basis
    Pi(x) = sum(j) coefs[i,j]x^j
    """
    result = 0
    for j in range(3):
      result += self.coefs[i,j] * x**j
    result *= self.Pol[i]
    return result

  
  def EvaluateMonomialDerivative(self, x, i):
    if i == 0:
      return 0.0
    elif i == 1:
      return self.coefs[1,1]
    elif i == 2:
      return 2.0 * self.coefs[2,2] * x
    else:
      raise Exception("max order : 2") 

  
  def EvaluateMonomialSecondDerivative(self, x, i):
    if i < 2 :
      return 0
    elif i == 2:
      return 2.0 * self.coefs[2,2]
    else:
      raise Exception("max order : 2")

  def EvaluateCurvature(x):
    result = 0
    for i in range(3):
      result += self.EvaluateMonomialSecondDerivative(x,i) * self.Pol[i]
    return result

  



class LegendreMultivariate(Polynomial):

  """
  2D Legendre polynomials of total order 2
  uniCoefs defines the monomials used in each dimension for a given basis polynomial
  if 0, then no monomial, if 1, monomial of order 0, if 2, monomial of order 1, etc.
  """
  
  Pol           = None
  xUni          = None
  nCoefs        = None
  uniCoefs      = None

  def __init__(self, Pol, dim):
    """
    Pol is a list of the coefficients for the monomials. It defines completely the polynomials 
    """
    self.dim = dim
    self.xUni = []
    if dim == 2:
      self.nCoefs = 6
      self.uniCoefs = np.zeros((self.nCoefs,self.dim), dtype=np.intc)
      self.uniCoefs[0,0] = 1 
      self.uniCoefs[0,1] = 1
      self.uniCoefs[1,0] = 2 
      self.uniCoefs[1,1] = 0
      self.uniCoefs[2,0] = 0 
      self.uniCoefs[2,1] = 2
      self.uniCoefs[3,0] = 3 
      self.uniCoefs[3,1] = 0
      self.uniCoefs[4,0] = 2 
      self.uniCoefs[4,1] = 2
      self.uniCoefs[5,0] = 0 
      self.uniCoefs[5,1] = 3

    elif dim == 3:
      self.nCoefs = 10
      self.uniCoefs = np.zeros((self.nCoefs,self.dim), dtype=np.intc)
      self.uniCoefs[0,0] = 1 
      self.uniCoefs[0,1] = 1
      self.uniCoefs[0,2] = 1
      self.uniCoefs[1,0] = 2 
      self.uniCoefs[1,1] = 0
      self.uniCoefs[1,2] = 0
      self.uniCoefs[2,0] = 0 
      self.uniCoefs[2,1] = 2
      self.uniCoefs[2,2] = 0
      self.uniCoefs[3,0] = 0 
      self.uniCoefs[3,1] = 0
      self.uniCoefs[3,2] = 2
      self.uniCoefs[4,0] = 3 
      self.uniCoefs[4,1] = 0
      self.uniCoefs[4,2] = 0
      self.uniCoefs[5,0] = 2 
      self.uniCoefs[5,1] = 2
      self.uniCoefs[5,2] = 0
      self.uniCoefs[6,0] = 2 
      self.uniCoefs[6,1] = 0
      self.uniCoefs[6,2] = 2
      self.uniCoefs[7,0] = 0 
      self.uniCoefs[7,1] = 3
      self.uniCoefs[7,2] = 0
      self.uniCoefs[8,0] = 0 
      self.uniCoefs[8,1] = 2
      self.uniCoefs[8,2] = 2
      self.uniCoefs[9,0] = 0 
      self.uniCoefs[9,1] = 0
      self.uniCoefs[9,2] = 3

    if len(Pol) == self.nCoefs:
      self.Pol = Pol
    else:
      raise Exception(str(self.nCoefs) + " coeffs are needed to define a " + str(dim) + "D polynomial of order 2")
    for i in range(self.dim): 
      self.xUni.append(LegendreUnivariate([1,1,1]))
 
 
  def Evaluate(self, x):
    result = 0
    for i in range(self.nCoefs):
      result += self.EvaluateMonomial(x, i)
    return result
 

  def EvaluateMonomial(self, x, n):
    """
    x : evaluation point
    n : index of the monomial in the basis
    """
    result = 0
    result = self.Pol[n]
    for j in range(self.dim):
      if self.uniCoefs[n,j] != 0:
        result *= self.xUni[j].EvaluateMonomial(x[j], self.uniCoefs[n,j]-1)
    return result


  def EvaluateDerivative(self, x, axis):
    result = 0
    for i in range(self.nCoefs):
      if self.uniCoefs[i, axis] != 0:
        temp = self.Pol[i]
        for j in range(self.dim):
          if j == axis:
            if self.uniCoefs[i,j] == 0:
              temp = 0
            else:
              temp *= self.xUni[j].EvaluateMonomialDerivative(x[j], self.uniCoefs[i,j]-1)
          else:
            if self.uniCoefs[i,j] == 0:
              pass # Nothing to do here
            else:
              temp *= self.xUni[j].EvaluateMonomial(x[j], self.uniCoefs[i,j]-1)
        result += temp
    return result


  def EvaluateSecondDerivative(self, x, axis1, axis2):
    result = 0
    if axis1 != axis2:
      for i in range(self.nCoefs):
        temp = self.Pol[i]
        for j in range(self.dim):
          if j == axis1:
            if self.uniCoefs[i,j] == 0:
              temp = 0
            else:
              temp *= self.xUni[j].EvaluateMonomialDerivative(x[j], self.uniCoefs[i,j]-1)
          if j == axis2:
            if self.uniCoefs[i,j] == 0:
              temp = 0
            else:
              temp *= self.xUni[j].EvaluateMonomialDerivative(x[j], self.uniCoefs[i,j]-1)
          else:
            if self.uniCoefs[i,j] == 0:
              pass # Nothing to do here
            else:
              temp *= self.xUni[j].EvaluateMonomial(x[j], self.uniCoefs[i,j]-1)
        result += temp
    else:
      for i in range(self.nCoefs):
        if self.uniCoefs[i, axis1] != 0:
          temp = self.Pol[i]
          for j in range(self.dim):
            if j == axis1:
              if self.uniCoefs[i,j] == 0:
                temp = 0
              else:
                temp *= self.xUni[j].EvaluateMonomialSecondDerivative(x[j], self.uniCoefs[i,j]-1)
            else:
              if self.uniCoefs[i,j] == 0:
                pass # Nothing to do here
              else:
                temp *= self.xUni[j].EvaluateMonomial(x[j], self.uniCoefs[i,j]-1)
          result += temp
    return result


  def EvaluateGradient(self,x):
    fx    = self.EvaluateDerivative(x,0)
    fy    = self.EvaluateDerivative(x,1)
    if self.dim == 3:
      fz  = self.EvaluateDerivative(x,2)
      return [fx,fy,fz]
    else:
      return [fx,fy]


  def EvaluateCurvature(self,x):
    one = 1.0
    two = 2.0
    epsilon = 1e-15
    fx    = self.EvaluateDerivative(x,0)
    fy    = self.EvaluateDerivative(x,1)
    if self.dim == 3:
      fz  = self.EvaluateDerivative(x,2)
    fxx   = self.EvaluateSecondDerivative(x,0,0)
    fyy   = self.EvaluateSecondDerivative(x,1,1)
    fxy   = self.EvaluateSecondDerivative(x,0,1)
    if self.dim == 3:
      fzz = self.EvaluateSecondDerivative(x,2,2)
      fzx = self.EvaluateSecondDerivative(x,2,0)
      fyz = self.EvaluateSecondDerivative(x,1,2)

    if self.dim == 2:
      norm  = math.sqrt( fx**2 + fy**2)
      nx    = fx/norm
      ny    = fy/norm
      kappa = -((one-nx**2)*fxx+(one-ny**2)*fyy- \
                                   two*nx*ny*fxy)/max(norm,epsilon)
    elif self.dim ==3:
      norm  = math.sqrt( fx**2 + fy**2 + fz**2 )
      nx    = fx/norm
      ny    = fy/norm
      nz    = fz/norm
      kappa = -((one-nx**2)*fxx+(one-ny**2)*fyy+(one-nz**2)*fzz- \
                                   two*nx*ny*fxy-two*nx*nz*fzx-two*ny*nz*fyz)/ \
                                   max(norm,epsilon)
    else:
      raise Exception("Curvature error, this should not happen")
    return kappa
