import os
import numpy as np
import matplotlib
#matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc
from matplotlib.colors import ListedColormap
import AVBP_Tools.AxiFrom3D.InterpolateFromAVBP as ci
from AVBP_Tools.AxiFrom3D.DetectAutoIgnition import DetectAutoIgnition as det
import common_IO as cio
from scipy.interpolate import interp1d
import Utils.PrettyPlots as pp

def GetField(workDir,variable, interpRelPath):
  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'av_sol'

  # Keyname of variable to integrate (as specified in solution h5-file)
  VARKEY = 'Average/' + variable

  ##########################################################################################

  # Read data
  interp_base = cio.load_interp_base(workDir + '/' + interpRelPath)
  wts = interp_base['wts']
  vtx = interp_base['vtx']
  has_voln = interp_base['header'][1][6]
  
  # Read bounding box of interpolation mesh from header
  xlim = np.array(interp_base['header'][1][3])
  rlim = np.array(interp_base['header'][1][4])
  tlim = np.array(interp_base['header'][1][5])
  
  
  # Get mesh resolution from header
  resx = interp_base['header'][1][0]
  resr = interp_base['header'][1][1]
  rest = interp_base['header'][1][2]
  
  
  # Generate interpolation mesh-vectors
  print "Generate interpolation mesh axes"
  xq = np.linspace(xlim[0],xlim[1],resx)
  rq = np.linspace(rlim[0],rlim[1],resr)
  tq = np.linspace(tlim[0],tlim[1],rest+1)[:-1] # We take every position except 360 degrees which is a duplicate of 0 degrees

  #print "Interpolation mesh contains %i cells" % (resx*resr*rest)
  
  # Loop through all solution files and apply interpolation on coarse grid





  counter = 0
  for solutionfile in os.listdir(workDir):
    if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
      assert counter==0,"Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

      hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),VARKEY)
      print "Interpolating"
      vq = ci.interpolate(hr,vtx,wts)

      print "axi mean"
      xQtemp = []
      for it in range(rest):
        imin = it * resx*resr
        imax = (it + 1) * resx*resr
        xQtemp.append(vq[imin:imax])
      xQtemp = np.array(xQtemp)
      mean = np.zeros(resx*resr)


      for ip in range(resx*resr):
        mean[ip] = np.sum(xQtemp[:,ip])/rest

  
   
  return mean, resx, resr, xlim, rlim
