import math
from Utils.I_O import ErrorMsg

fileName = 'VarTransformation.py'

# Direct transformations

def NormalToLogNormal(value, xTransformationParam):
  nominalValue  = xTransformationParam[0]
  UF            = xTransformationParam[1]
  factor        = math.exp(math.log(UF) * value / 3)
  result        = nominalValue * factor
  return result

def NoTransform(value, xTransformationParam):
  return value

def Dirac(value, xTransformationParam):
  return xTransformationParam[0]

def StretchedUniform(value, xTransformationParam):
  inf     = xTransformationParam[0]
  sup     = xTransformationParam[1]
  delta   = sup - inf
  result  = value * delta + inf
  return result

def UniformToLogUniform(value, xTransformationParam):
  inf     = math.log(xTransformationParam[0])
  sup     = math.log(xTransformationParam[1])
  delta   = sup - inf
  result  = math.exp(value * delta + inf)
  return result


# Inverse transformations

def NormalToLogNormalInverse(value, xTransformationParam):
  nominalValue  = xTransformationParam[0]
  UF            = xTransformationParam[1]
  result        = 3 * math.log(value / nominalValue - UF)
  return result

def NoTransformInverse(value, xTransformationParam):
  return value

def DiracInverse(value, xTransformationParam):
  ErrorMsg('This transformation is not invertible', fileName)
  return None

def StretchedUniformInverse(value, xTransformationParam):
  inf     = xTransformationParam[0]
  sup     = xTransformationParam[1]
  delta   = sup - inf
  result  = (value - inf) / delta
  return result

def UniformToLogUniformInverse(value, xTransformationParam):
  inf     = math.log(xTransformationParam[0])
  sup     = math.log(xTransformationParam[1])
  delta   = sup - inf
  result  = (math.log(value) + inf) / delta
  return result
