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

def PlotField(workDir, D, variable, interpRelPath, solRelPath=None, xIso=[], bottomThresh=None):
  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'av_sol'

  # Save images to separate folder
  IMG_OUT = workDir

  # Subdivisions of axis
  xstep = 10
  rstep = 2

  # Keyname of variable to integrate (as specified in solution h5-file)
  VARKEY = 'Average/' + variable

  if variable=='T':
    cLab = r'Temperature [K]'
  elif variable == 'Yzv':
    cLab = r'Zv [-]'
  elif variable == 'z_out':
    cLab = r'Z [-]'
  elif variable=='c_out':
    #cLab = r'C [-]'
    cLab = r'Normalized Progress Variable [-]'
  elif variable == 'WYc_out':
    cLab = r'$\dot{\omega}_{Y_C}$ [s$^{-1}$]'
  else:
    print('Please check that that the y label corresponds to what you expect')
    cLab = r'Y$_{' + variable + '}$ [-]'

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
      #assert counter==0,"Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

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

      Q = np.reshape(mean,(resr,resx))

      Q[np.isnan(Q)] = 0


      # Plot
      print "Prepare plot and save figure"
      fig, ax = plt.subplots(1, 1, figsize=(7, 7))
      # font = {'family' : 'sans',
      #    'serif': 'helvet',
      #    'weight' : 'bold',
      #    'size'   : 16}
      #
      # plt.rc('font', **font)
      # plt.rc('text', usetex=True)
      #cmap = 'hot'

      #Custom colormap to match Cabra 2002 for OH

      def forceAspect(ax,aspect):
        im = ax.get_images()
        extent =  im[0].get_extent()
        ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

      if variable == 'OH':
        norm = matplotlib.colors.Normalize(vmin=0, vmax=0.001) # for OH
        cmap = ListedColormap([(0.15,0.18,0.35),(0.15,0.2,0.38),(0.2,0.36,0.58),(0.23,0.42,0.64),(0.25,0.51,0.65),(0.21,0.4,0.24),(0.25,0.49,0.28),(0.54,0.71,0.33),(0.97,0.91,0.37),(0.73,0.33,0.23),(0.67,0.22,0.2),(0.72,0.18,0.2)])
      elif variable == 'WYc_out':
        norm = matplotlib.colors.Normalize(vmin=0, vmax=30)
        #cmap = 'Blues'
        cmap = 'hot_r'
        cmap = 'jet'
        #cmap = 'Reds'
      elif variable == 'T':
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1500) # for T
        cmap = 'hot'
      elif variable == 'YZv':
        norm = matplotlib.colors.Normalize(vmin=0, vmax=0.01) # for zv
        cmap = 'plasma'
      else:
        norm = matplotlib.colors.Normalize(vmin=0, vmax=Q.max())
        cmap='viridis'
      #norm = matplotlib.colors.Normalize(vmin=0, vmax=1500)
      #norm = matplotlib.colors.Normalize(vmin=0, vmax=0.004)
      ax = plt.gca()
      img = ax.imshow(np.rot90(Q, 1), interpolation='bilinear', cmap=cmap, norm=norm,
                      extent=[rq.min(), rq.max(), xq.min(), xq.max()])
      forceAspect(ax, aspect=0.25)
      #plt.gca().set_aspect('equal', adjustable='box')
      #plt.gca().set_aspect(4.0, adjustable='box')
      #fig.set_figwidth(4)
      #fig.set_figheight(8)



      if len(xIso)>0:
        myContour = plt.contour(rq, xq, np.flipud(np.rot90(Q, 1)), xIso, colors='w', linestyles='--', linewidths=1.5)

      cbar = plt.colorbar(img)
      cbar.set_label(cLab, size=30)
      cbar.ax.tick_params(labelsize=30)

      xticks = np.arange(rlim[0], rlim[1], rstep * D)
      xtickslabels = np.round(xticks / D)
      plt.xticks(xticks, xtickslabels, fontsize=30)

      yticks = np.arange(xlim[0], xlim[1], xstep * D)
      ytickslabels = np.round(yticks / D)
      plt.yticks(yticks, ytickslabels, fontsize=30)

      plt.xlabel('r/D', fontsize=30)
      plt.ylabel('x/D', fontsize=30)

      pp.CorrectAxSize(ax, 30)

      if bottomThresh is not None:
        altqoi = det(Q,bottomThresh, xq,rq)/D
        print(altqoi)
      plt.savefig(workDir + '/Field_' + variable + '.pdf', bbox_inches="tight",pad_inches = 0)
      counter +=1

    else:
      continue
  
  
  
