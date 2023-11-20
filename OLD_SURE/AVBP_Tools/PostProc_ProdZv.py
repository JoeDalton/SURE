import os
import numpy as np
import matplotlib
#matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc
import AVBP_Tools.AxiFrom3D.InterpolateFromAVBP as ci
from AVBP_Tools.AxiFrom3D.DetectAutoIgnition import DetectAutoIgnition as det
import AVBP_Tools.AxiFrom3D.common_IO as cio
from scipy.interpolate import interp1d
#import Utils.PrettyPlots as pp

def PostProcProdZv(workDir, D, interpRelPath, solRelPath=None):
  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'current'

  # Save images to separate folder
  IMG_OUT = workDir

  # Subdivisions of axis
  xstep = 10
  ystep = 10


  variable = 'Z'

  # Keyname of variable to integrate (as specified in solution h5-file)
  VARKEY = 'FictiveSpecies/' + variable

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
  zq = np.linspace(zlim[0],zlim[1],resz) # We take every position except 360 degrees which is a duplicate of 0 degrees

  #print "Interpolation mesh contains %i cells" % (resx*resr*rest)
  dx = xq[1] - xq[0]
  dy = yq[1] - yq[0]
  dz = zq[1] - zq[0]
   
  # Loop through all solution files and apply interpolation on coarse grid
  counter = 0
  for solutionfile in os.listdir(workDir):
    if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
      assert counter==0,"Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

      print "Retrieving fields"
      Z       = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),'FictiveSpecies/Z')
      mu      = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),'Additionals/vis_turb')
      rho     = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),'GaseousPhase/rho')
      prod    = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),'Additionals/Prod_zv_out')
      print "Interpolating Z"
      Zq        = ci.interpolate(Z,vtx,wts)
      Zinterp   = np.reshape(Zq, (resy, resx, resz))
      print "Interpolating mu"
      muq       = ci.interpolate(mu,vtx,wts)
      muinterp  = np.reshape(muq, (resy, resx, resz))
      print "Interpolating rho"
      rhoq      = ci.interpolate(rho,vtx,wts)
      rhointerp = np.reshape(rhoq, (resy, resx, resz))
      print "Interpolating Prod_zv_out"
      pq        = ci.interpolate(prod,vtx,wts)
      pinterp   = np.reshape(pq, (resy, resx, resz))
      

      """
      Operators
      """

      def ComputeGrad(Field, resx, resy, resz, dx, dy, dz):
        gradX  = np.zeros_like(Field)
        gradY  = np.zeros_like(Field)
        gradZ  = np.zeros_like(Field)
        #Computing dZ/dX
        gradX[:,0,:]       = (Field[:,1,:]      - Field[:,0,:])        / (dx)
        gradX[:,resx-1,:]  = (Field[:,resx-1,:] - Field[:,resx-2,:])    / (dx)
        for i in range(1, resx-1):
          gradX[:,i,:]     = (Field[:,i+1,:]    - Field[:,i-1,:])      / (2.0*dx)
        
        #Computing dZ/dY
        gradY[0,:,:]       = (Field[1,:,:]      - Field[0,:,:])        / (dy)
        gradY[resy-1,:,:]  = (Field[resy-1,:,:] - Field[resy-2,:,:])   / (dy)
        for i in range(1, resy-1):
          gradY[i,:,:]     = (Field[i+1,:,:]    - Field[i-1,:,:])      / (2.0*dy)
        
        #Computing dZ/dZ
        gradZ[:,:,0]       = (Field[:,:,1]      - Field[:,:,0])        / (dz)
        gradZ[:,:,resz-1]  = (Field[:,:,resz-1] - Field[:,:,resz-2])   / (dz)
        for i in range(1, resz-1):
          gradZ[:,:,i]     = (Field[:,:,i+1]    - Field[:,:,i-1])      / (2.0*dz)
        
        return gradX, gradY, gradZ


      def ComputeDiv(VectX, VectY, VectZ, resx, resy, resz, dx, dy, dz):
        divX  = np.zeros_like(VectX)
        divY  = np.zeros_like(VectX)
        divZ  = np.zeros_like(VectX)
        div   = np.zeros_like(VectX)
        #Computing dZx/dX
        divX[:,0,:]       = (VectX[:,1,:]      - VectX[:,0,:])        / (dx)
        divX[:,resx-1,:]  = (VectX[:,resx-1,:] - VectX[:,resx-2,:])    / (dx)
        for i in range(1, resx-1):
          divX[:,i,:]     = (VectX[:,i+1,:]    - VectX[:,i-1,:])      / (2.0*dx)
        
        #Computing dZy/dY
        divY[0,:,:]       = (VectY[1,:,:]      - VectY[0,:,:])        / (dy)
        divY[resy-1,:,:]  = (VectY[resy-1,:,:] - VectY[resy-2,:,:])   / (dy)
        for i in range(1, resy-1):
          divY[i,:,:]     = (VectY[i+1,:,:]    - VectY[i-1,:,:])      / (2.0*dy)
        
        #Computing dZz/dZ
        divZ[:,:,0]       = (VectZ[:,:,1]      - VectZ[:,:,0])        / (dz)
        divZ[:,:,resz-1]  = (VectZ[:,:,resz-1] - VectZ[:,:,resz-2])   / (dz)
        for i in range(1, resz-1):
          divZ[:,:,i]     = (VectZ[:,:,i+1]    - VectZ[:,:,i-1])      / (2.0*dz)
        
        div[:,:,:] = divX[:,:,:] + divY[:,:,:] + divZ[:,:,:]
        return div


      """
      "Naive" |grad(Z)|^2
      """
      gradZX  = np.zeros_like(Zinterp)
      gradZY  = np.zeros_like(Zinterp)
      gradZZ  = np.zeros_like(Zinterp)
      grad2Zn = np.zeros_like(Zinterp)

      gradZX, gradZY, gradZZ = ComputeGrad(Zinterp, resx, resy, resz, dx, dy, dz)

      grad2Zn[:,:,:]      = gradZX[:,:,:] * gradZX[:,:,:] + gradZY[:,:,:] * gradZY[:,:,:] + gradZZ[:,:,:] * gradZZ[:,:,:]

      """
      "Smart" |grad(Z)|^2
      """
      Z2                  = np.zeros_like(Zinterp)
      Z2[:,:,:]           = Zinterp[:,:,:] * Zinterp[:,:,:]
      gradZ2X, gradZ2Y, gradZ2Z = ComputeGrad(Z2, resx, resy, resz, dx, dy, dz)
      DivGradZ2 = ComputeDiv(gradZ2X, gradZ2Y, gradZ2Z, resx, resy, resz, dx, dy, dz)
      DivGradZ  = ComputeDiv(gradZX, gradZY, gradZZ, resx, resy, resz, dx, dy, dz)

      grad2Zs = np.zeros_like(Zinterp)
      grad2Zs[:,:,:] = 0.5 * (DivGradZ2[:,:,:] - 2.0 * Zinterp[:,:,:] * DivGradZ[:,:,:])




      """
      PostProcessing Prod_zv
      """
      Schmidt = 0.6

      rhoD_t  = np.zeros_like(Zinterp)
      ppp     = np.zeros_like(Zinterp)

      rhoD_t[:,:,:] = muinterp[:,:,:] / Schmidt
      #ppp[:,:,:]    = 2.0 * rhoD_t[:,:,:] * grad2Zn[:,:,:] / rhointerp[:,:,:]
      ppp[:,:,:]    = 2.0 * rhoD_t[:,:,:] * grad2Zs[:,:,:] / rhointerp[:,:,:]
      

      #Q = Zinterp[:,:,int(abs(resz/2))]
      #Q = grad2Zn[:,:,int(abs(resz/2))]
      #Q = grad2Zs[:,:,int(abs(resz/2))]

      Q1 = pinterp[:,:,int(abs(resz/2))]
      Q2 = ppp[:,:,int(abs(resz/2))]


      # Plot

      
      cLab = r'Prod of zv'


      print "Prepare plot and save figure"
      #fig, ax = plt.subplots(1, 1, figsize=(12, 12))
      fig = plt.figure(1, figsize=(14, 10))
      """
      Plotting initial field
      """
      ax1 = fig.add_subplot(121)
      #cmap = 'hot'
      cmap = 'Blues'
      #norm = matplotlib.colors.Normalize(vmin=0, vmax=1500)
      #norm1 = matplotlib.colors.Normalize(vmin=0, vmax=Q1.max())
      #norm2 = matplotlib.colors.Normalize(vmin=0, vmax=Q2.max())
      norm1 = matplotlib.colors.Normalize(vmin=0, vmax=4000)
      norm2 = matplotlib.colors.Normalize(vmin=0, vmax=300)
      #norm = matplotlib.colors.Normalize(vmin=0, vmax=1.0e6)
      #norm = matplotlib.colors.Normalize(vmin=0, vmax=1.0e16)
      #ax = plt.gca()
      img1 = ax1.imshow(np.rot90(Q1, 1), interpolation='bilinear', cmap=cmap, norm=norm1,
                      extent=[yq.min(), yq.max(), xq.min(), xq.max()])

      cbar1 = plt.colorbar(img1)
      cbar1.set_label(cLab, size=16)

      ax1.set_title('From AVBP')
      """
      Plotting post proc
      """
      ax2 = fig.add_subplot(122)
      # font = {'family' : 'sans',
      #    'serif': 'helvet',
      #    'weight' : 'bold',
      #    'size'   : 16}
      #
      # plt.rc('font', **font)
      # plt.rc('text', usetex=True)
      #cmap = 'hot'
      cmap = 'Blues'
      img2 = ax2.imshow(np.rot90(Q2, 1), interpolation='bilinear', cmap=cmap, norm=norm2,
                      extent=[yq.min(), yq.max(), xq.min(), xq.max()])
      cbar2 = plt.colorbar(img2)
      cbar2.set_label(cLab, size=16)
      ax2.set_title('Post-proc')
      plt.gca().set_aspect('equal', adjustable='box')

      #cbar.ax1.tick_params(labelsize=12)

      xticks = [-10.0*D, 0.0, 10.0*D]
      xtickslabels = [-10.0, 0.0, 10.0]
      #xtickslabels = np.round(xticks / D)
      #plt.xticks(xticks, xtickslabels, fontsize=12)
      ax1.set_xticks(xticks)
      ax1.set_xticklabels(xtickslabels)
      ax2.set_xticks(xticks)
      ax2.set_xticklabels(xtickslabels)

      yticks = np.arange(xlim[0], xlim[1], xstep * D)
      ytickslabels = np.round(yticks / D)
      #plt.yticks(yticks, ytickslabels, fontsize=12)
      ax1.set_yticks(yticks)
      ax1.set_yticklabels(ytickslabels)
      ax2.set_yticks(yticks)
      ax2.set_yticklabels(ytickslabels)

      plt.xlabel('y/D', fontsize=16)
      plt.ylabel('x/D', fontsize=16)

      #plt.show()
      plt.savefig('prod_smart.pdf', bbox_inches="tight",pad_inches = 0)
      counter +=1

    else:
      continue
  
  
  
