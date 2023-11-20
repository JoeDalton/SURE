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

def GetInstantaneous(workDir, variable, interpRelPath):

  # interpRelPath must designate a cartesian mesh with 3 points in the direction z (slice at z = 0)

  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'inst'

  # Keyname of variable to integrate (as specified in solution h5-file)
  VARKEY = 'Additionals/' + variable

  ##########################################################################################

  # Read data
  interp_base = cio.load_interp_base(workDir + '/' + interpRelPath)
  wts = interp_base['wts']
  vtx = interp_base['vtx']
  has_voln = interp_base['header'][1][6]
  
  # Read bounding box of interpolation mesh from header
  xlim = np.array(interp_base['header'][1][3])
  ylim = np.array(interp_base['header'][1][4])
  zlim = np.array(interp_base['header'][1][5])
  
  
  # Get mesh resolution from header
  resx = interp_base['header'][1][0]
  resy = interp_base['header'][1][1]
  resz = interp_base['header'][1][2]
  
  
  # Generate interpolation mesh-vectors
  print "Generate interpolation mesh axes"
  xq = np.linspace(xlim[0],xlim[1],resx)
  yq = np.linspace(ylim[0],ylim[1],resy)
  zq = np.linspace(zlim[0],zlim[1],resz)

  #print "Interpolation mesh contains %i cells" % (resx*resr*rest)
  
  # Loop through all solution files and apply interpolation on coarse grid





  counter = 0
  solutionList = sorted(os.listdir(workDir))
  for solutionfile in solutionList:
    if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
      assert counter==0,"Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

      hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),VARKEY)
      print "Interpolating"
      vq = ci.interpolate(hr,vtx,wts)
      inst=vq

      #print "axi mean"
      #xQtemp = []
      #for it in range(rest):
      #  imin = it * resx*resr
      #  imax = (it + 1) * resx*resr
      #  xQtemp.append(vq[imin:imax])
      #xQtemp = np.array(xQtemp)
      #mean = np.zeros(resx*resr)


      #for ip in range(resx*resr):
      #  mean[ip] = np.sum(xQtemp[:,ip])/rest

  
   
  return inst, resx, resy, resz, xlim, ylim, zlim




def GetInstantaneousFrom2DSlice(workDir, variable, interpRelPath):

  # interpRelPath must designate a cartesian mesh with 3 points in the direction z (slice at z = 0)

  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'mycut_cut'

  # Keyname of variable to integrate (as specified in solution h5-file)
  VARKEY = 'Additionals/' + variable

  ##########################################################################################

  # Read data
  interp_base = cio.load_interp_base(workDir + '/' + interpRelPath)
  wts = interp_base['wts']
  vtx = interp_base['vtx']
  
  # Read bounding box of interpolation mesh from header
  xlim = np.array(interp_base['header'][1][2])
  ylim = np.array(interp_base['header'][1][3])
  
  
  # Get mesh resolution from header
  resx = interp_base['header'][1][0]
  resy = interp_base['header'][1][1]
  
  
  # Generate interpolation mesh-vectors
  print "Generate interpolation mesh axes"
  xq = np.linspace(xlim[0],xlim[1],resx)
  yq = np.linspace(ylim[0],ylim[1],resy)

  #print "Interpolation mesh contains %i cells" % (resx*resr*rest)
  
  # Loop through all solution files and apply interpolation on coarse grid



  xInst = []

  solutionList = sorted(os.listdir(workDir))
  for solutionfile in solutionList:
    if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
      #print(solutionfile)
      hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),VARKEY)
      print "Interpolating"
      vq = ci.interpolate(hr,vtx,wts)
      xInst.append(vq)

      #print "axi mean"
      #xQtemp = []
      #for it in range(rest):
      #  imin = it * resx*resr
      #  imax = (it + 1) * resx*resr
      #  xQtemp.append(vq[imin:imax])
      #xQtemp = np.array(xQtemp)
      #mean = np.zeros(resx*resr)


      #for ip in range(resx*resr):
      #  mean[ip] = np.sum(xQtemp[:,ip])/rest

  
   
  return xInst, resx, resy, xlim, ylim
