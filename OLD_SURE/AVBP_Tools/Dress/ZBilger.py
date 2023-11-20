import os
import numpy as np
import AVBP_Tools.AxiFrom3D.common_IO as cio





def DressZBilgerCabra(solPath):
  # Define molar mass of atoms
  MH          = 1.0
  MO          = 16.0

  # Define mass fractions of atoms in the jet and coflow
  YHfu        = 0.0239923
  YHcof       = 0.007171833
  YOfu        = 0.0
  YOcof       = 0.2282443

  # Read mass fraction of species from the solution
  YH2         = cio.read_avbp_solution(solPath,"H2")
  YO2         = cio.read_avbp_solution(solPath,"O2")
  YH2O        = cio.read_avbp_solution(solPath,"H2O")
  YOH         = cio.read_avbp_solution(solPath,"OH")
  YH          = cio.read_avbp_solution(solPath,"H")
  YHO2        = cio.read_avbp_solution(solPath,"HO2")
  YO          = cio.read_avbp_solution(solPath,"O")

  # Define mass fraction of atoms in molecules
  YH_H2       = 1.0 
  YH_O2       = 0.0
  YH_H2O      = 1.0/9.0
  YH_OH       = 1.0/17.0
  YH_H        = 1.0
  YH_HO2      = 1.0/33.0
  YH_O        = 0.0

  YO_H2       = 0.0
  YO_O2       = 1.0
  YO_H2O      = 8.0/9.0
  YO_OH       = 16.0/17.0
  YO_H        = 0.0
  YO_HO2      = 32.0/33.0
  YO_O        = 1.0

  # Initialize fields
  ZBilger     = np.zeros_like(YH2)
  num         = np.zeros_like(YH2)
  YOat        = np.zeros_like(YH2)
  YHat        = np.zeros_like(YH2)
 
  # Compute Fields
  YHat[:]     = YH_H2 * YH2[:] + YH_O2 * YO2[:] + YH_H2O * YH2O[:] + YH_OH * YOH[:] + YH_H * YH[:] + YH_HO2 * YHO2[:] + YH_O * YO[:]
  YOat[:]     = YO_H2 * YH2[:] + YO_O2 * YO2[:] + YO_H2O * YH2O[:] + YO_OH * YOH[:] + YO_H * YH[:] + YO_HO2 * YHO2[:] + YO_O * YO[:]

  num[:]      = (YHat[:]  - YHcof) / (2 * MH) - (YOat[:]  - YOcof) / MO
  denom       = (YHfu     - YHcof) / (2 * MH) - (YOfu     - YOcof) / MO
  print num.min()
  
  ZBilger[:]  = num[:] / denom

  # Write ZBilger field
  cio.write_avbp_solution(solPath, "Average/ZBilger", ZBilger)






