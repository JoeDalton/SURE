#!/usr/bin/env python2.7

from __future__ import print_function
import h5py
import time
import Link
import numpy                                      as np
import openturns                                  as ot
import matplotlib.pyplot                          as plt
from SurrogateModelling.PCE             import *
from SurrogateModelling.OneDSpline      import *
from DimensionReduction.ActiveDirection import *
from DimensionReduction.Preprocessing import *
import Utils.PrettyPlots                          as pp

t0 = time.time()

nSample = 5000
#nSample = 20


counter = 0

for T in [913.5]:
  #############################################################
  # STEP -1: Retrieving the sampling data
  #############################################################
  tFolderName	= './testData/T_' + str(T)
  print(tFolderName)
  resultFile = tFolderName + '/result.h5'
  result  = h5py.File(resultFile, 'r')
  xEval = np.log(result["Data/QoI"][:nSample])
  xSample = result["Data/NormalInputVars"][:nSample,:]
  result.close()
  
  ##Uncertain dimension and size of the dataset
  #dim = len(points[0,:])
  #dataSize = len(evaluations)
  
  
  #############################################################
  # STEP 0: Creating the Active Direction object
  #############################################################
  objAD = ActiveDirection(ID='tata', totalDim=xSample.shape[1], verbosity=1)
  objAD.Build(xSample, xEval)
  objAD.Save()

  #############################################################
  # STEP 1: Constructing the list of reduced samples
  #############################################################
  xRedSample = objAD.ComputeReducedFromTotal(xSample)

  #############################################################
  # STEP 2: Centering and scaling the data
  #############################################################
  myPreProc = CenterScaling("toto", xEval, "RANGE")
  xNormEval = myPreProc.ComputePro(xEval)

  #############################################################
  # STEP 3: Building a spline surrogate from the reduced var
  #############################################################
  objSpl = OneDSpline(ID='toto', verbosity=1)
  objSpl.Build(samples=xRedSample, evaluations=xEval, trainingRatio=0.8)
  xReEval = np.zeros_like(xRedSample)
  for i in range(len(xRedSample)):
    xReEval[i] = objSpl.Evaluate(xRedSample[i])
  xNormReEval = myPreProc.ComputePro(xReEval)


  #############################################################
  # STEP 3bis: CPPN Surrogates
  #############################################################
  def CPPN_1(Inputs):
    # CPPN GuyId:1304100461040 EspeceId:1504100458604
    assert len(Inputs) == 1, 'Nbr d entree incorrect'
    Outputs = [None] * 1

    # NodeId:1104100461681 - Type:INPUT - Aggregation:INPUT - Transfert:IDENT
    N_1104100461681_OUT = (Inputs[0])

    # NodeId:1104100461683 - Type:HIDDEN - Aggregation:ADD - Transfert:IDENT
    N_1104100461683_Aggreg = 0.792450
    N_1104100461683_Aggreg += N_1104100461681_OUT * (-3.202824)
    N_1104100461683_OUT = (N_1104100461683_Aggreg)

    # NodeId:1104100461682 - Type:OUTPUT - Aggregation:ADD - Transfert:TANH
    N_1104100461682_Aggreg = -0.731272
    N_1104100461682_Aggreg += N_1104100461681_OUT * (2.792103)
    N_1104100461682_Aggreg += N_1104100461683_OUT * (0.805912)
    N_1104100461682_OUT = math.tanh(N_1104100461682_Aggreg)
    Outputs[0] = N_1104100461682_OUT

    return Outputs
  def SOFT_RELU (x):
      if x > 20.0:return x
      else:       return math.log1p ( 1 + math.exp(x))

  def CPPN_2(Inputs):
    # CPPN GuyId:1304100479563 EspeceId:1504100458604
    assert len(Inputs) == 1, 'Nbr d entree incorrect'
    Outputs = [None] * 1

    # NodeId:1104100504218 - Type:INPUT - Aggregation:INPUT - Transfert:IDENT
    N_1104100504218_OUT = (Inputs[0])

    # NodeId:1104100528848 - Type:HIDDEN - Aggregation:ADD - Transfert:SOFT_RELU
    N_1104100528848_Aggreg = -0.488692
    N_1104100528848_Aggreg += N_1104100504218_OUT * (1.000000)
    N_1104100528848_OUT = SOFT_RELU(N_1104100528848_Aggreg)

    # NodeId:1104100504220 - Type:HIDDEN - Aggregation:ADD - Transfert:IDENT
    N_1104100504220_Aggreg = 0.863555
    N_1104100504220_Aggreg += N_1104100504218_OUT * (-3.027425)
    N_1104100504220_OUT = (N_1104100504220_Aggreg)

    # NodeId:1104100513676 - Type:HIDDEN - Aggregation:ADD - Transfert:TANH
    N_1104100513676_Aggreg = -0.445923
    N_1104100513676_Aggreg += N_1104100504218_OUT * (0.555954)
    N_1104100513676_Aggreg += N_1104100528848_OUT * (0.055595)
    N_1104100513676_Aggreg += N_1104100504220_OUT * (-0.866020)
    N_1104100513676_OUT = math.tanh(N_1104100513676_Aggreg)

    # NodeId:1104100504221 - Type:HIDDEN - Aggregation:ADD - Transfert:TANH
    N_1104100504221_Aggreg = -0.073568
    N_1104100504221_Aggreg += N_1104100504218_OUT * (1.094776)
    N_1104100504221_Aggreg += N_1104100513676_OUT * (0.109478)
    N_1104100504221_OUT = math.tanh(N_1104100504221_Aggreg)

    # NodeId:1104100504219 - Type:OUTPUT - Aggregation:ADD - Transfert:TANH
    N_1104100504219_Aggreg = -0.680782
    N_1104100504219_Aggreg += N_1104100504218_OUT * (2.495457)
    N_1104100504219_Aggreg += N_1104100504220_OUT * (0.805912)
    N_1104100504219_Aggreg += N_1104100504221_OUT * (0.278163)
    N_1104100504219_OUT = math.tanh(N_1104100504219_Aggreg)
    Outputs[0] = N_1104100504219_OUT

    return Outputs

  xCPPN1Eval = np.zeros_like(xRedSample)
  xCPPN2Eval = np.zeros_like(xRedSample)
  for i in range(len(xRedSample)):
    xCPPN1Eval[i] = CPPN_1([xRedSample[i]])[0]
    xCPPN2Eval[i] = CPPN_2([xRedSample[i]])[0]

  #############################################################
  # STEP 4: Some plots
  #############################################################
  x2 = np.linspace(xRedSample.min(), xRedSample.max(), 1000)
  fig = plt.figure(1, figsize=(7,5))
  #pp.AdjustSubplotSpacings()
  ax3 = fig.add_subplot(111)
  ax3.scatter(xRedSample, xNormEval, color='k', alpha = 1.0, marker=',', s=1, label='Samples')
  ax3.scatter(xRedSample, xNormReEval, color='g', alpha = 1.0, marker=',', s=1, label='Spline')
  ax3.scatter(xRedSample, xCPPN1Eval, color='b', alpha = 1.0, marker=',', s=1, label='CPPN1')
  ax3.scatter(xRedSample, xCPPN2Eval, color='r', alpha = 1.0, marker=',', s=1, label='CPPN2')
  #ax3.plot(x2, objSpl.Evaluate(x2), color='g', label='Spline surrogate')
  ax3.set_xlabel("Reduced variable")
  ax3.set_ylabel("QoI")
  ax3.set_title("Model plot (Spline)")
  ax3.legend(loc='best')
  #pp.CorrectAxSize(ax3)

  plt.show()
  plt.clf()

  #############################################################
  # STEP 5: write to txt file
  #############################################################
  #np.savetxt("x.txt", xRedSample, newline=" ")
  #np.savetxt("y.txt", xNormEval, newline=" ")
