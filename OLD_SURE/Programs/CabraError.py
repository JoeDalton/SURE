
from Utils.I_O import prt, LoadCSV
import numpy as np
from scipy.interpolate import interp1d
import math


def ComputeCabraError(metaVect, resr, resx, rlim, xlim):
  # metaVect must be the concatenation of the vectorized fields of OH, O2, H2, H2O, T, z_out (in this order)  
  xVar    = ['OH', 'O2', 'H2', 'H2O', 'T', 'z_out']
  nVar    = len(xVar)
  xVect   = np.split(metaVect, nVar, axis=0)
  xq = np.linspace(xlim[0],xlim[1],resx)
  rq = np.linspace(rlim[0],rlim[1],resr)
  D =  0.00457

  csvDir  = 'ExpProfiles'
  xNX     = [1,9,11,14,26]

  sqError = 0.0

  for iVar in range(nVar):
    # Compute squared error for this field
    var           = xVar[iVar]
    if var == 'OH':
      vect = NormOH(xVect[iVar])
    elif var == 'O2':
      vect = NormO2(xVect[iVar])
    elif var == 'H2':
      vect = NormH2(xVect[iVar])
    elif var == 'H2O':
      vect = NormH2O(xVect[iVar])
    elif var == 'T':
      vect = NormT(xVect[iVar])
    elif var == 'z_out':
      vect = NormZ(xVect[iVar])
    else:
      print("Whoever coded this is worthless, sorry.")
      quit()


    field         = np.reshape(vect, (resr, resx))
    fieldSqError  = 0.0
    if var=="z_out":
      csvPrefx = csvDir + '/' + 'Z' + '_'
    else:
      csvPrefx = csvDir + '/' + var + '_'
    csvSufx = '.csv'

    # On the axis:
    csvName   = csvPrefx + 'Axial' + csvSufx
    xX,xYRef  = LoadCSV(csvName)
    if var == 'OH':
      xYRef = NormOH(xYRef)
    elif var == 'O2':
      xYRef = NormO2(xYRef)
    elif var == 'H2':
      xYRef = NormH2(xYRef)
    elif var == 'H2O':
      xYRef = NormH2O(xYRef)
    elif var == 'T':
      xYRef = NormT(xYRef)
    elif var == 'z_out':
      xYRef = NormZ(xYRef)
    else:
      print("You lose.")
      quit()
    interpAx  = interp1d(xq, field[0,:], kind='cubic')
    for j in range(len(xX)):
      yField        = interpAx(xX[j]*D)
      fieldSqError += (yField - xYRef[j]) * (yField - xYRef[j])

    # For radial data:
    for i in range(len(xNX)):
      csvName     = csvPrefx + 'z:d=' + str(xNX[i]) + csvSufx
      xX,xYRef    = LoadCSV(csvName)
      if var == 'OH':
        xYRef = NormOH(xYRef)
      elif var == 'O2':
        xYRef = NormO2(xYRef)
      elif var == 'H2':
        xYRef = NormH2(xYRef)
      elif var == 'H2O':
        xYRef = NormH2O(xYRef)
      elif var == 'T':
        xYRef = NormT(xYRef)
      elif var == 'z_out':
        xYRef = NormZ(xYRef)
      else:
        print("Too bad.")
        quit()
      profile     = []
      # build radial profile
      for k in range(resr):
        interpk = interp1d(xq, field[k,:], kind='cubic')
        profile.append(interpk(xNX[i]*D))
      # find radial values
      interpAx  = interp1d(rq, profile, kind='cubic')
      for j in range(len(xX)):
        yField        = interpAx(abs(xX[j])/1000) # Because some exp results are taken on the other side of the axis and abscissae were written in mm...
        fieldSqError += (yField - xYRef[j]) * (yField - xYRef[j])
      
    sqError += fieldSqError

  # Take square root of error and return
  error = math.sqrt(sqError)
  return error


#TODO: Define values for the norms

def NormOH(vect):
  norm = 1.0e-3
  return vect / norm

def NormO2(vect):
  norm = 0.17
  return vect / norm

def NormH2(vect):
  norm = 0.025
  return vect / norm

def NormH2O(vect):
  norm = 0.11
  return vect / norm

def NormT(vect):
  minVal = 300.0
  maxVal = 1500.0
  return (vect - minVal) / (maxVal - minVal)

def NormZ(vect):
  norm = 1.0
  return vect / norm
