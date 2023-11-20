import os
import numpy as np
import matplotlib
#matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import AVBP_Tools.AxiFrom3D.InterpolateFromAVBP as ci
import common_IO as cio
from scipy.interpolate import interp1d



def PlotProfile(workDir, D, variable, label, xAx, xNX, interpRelPath, solRelPath=None, axial=False):
  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'av_sol'

  # Save images to separate folder
  IMG_OUT = workDir
  
  assert len(xAx)==len(xNX), 'The number of inputed figs must be the same as the number of axial positions'
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



  if axial:
    ax = xAx[0]
    counter = 0
    for solutionfile in os.listdir(workDir):
      if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
        assert counter == 0, "Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

        hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir), solutionfile), VARKEY)
        print "Interpolating"
        vq = ci.interpolate(hr, vtx, wts)

        print "axi mean"
        xQtemp = []
        for it in range(rest):
          imin = it * resx * resr
          imax = (it + 1) * resx * resr
          xQtemp.append(vq[imin:imax])
        xQtemp = np.array(xQtemp)
        mean = np.zeros(resx * resr)

        for ip in range(resx * resr):
          mean[ip] = np.sum(xQtemp[:, ip]) / rest

        Q = np.reshape(mean, (resr, resx))

        Q[np.isnan(Q)] = 0

        ax.plot(xq, Q[0,:], label=label)
        ax.legend(fontsize=12)
        counter += 1
      else:
        continue

  else:
    for i in range(len(xNX)):
      normAxialDistance = xNX[i]
      #fig, ax = plt.subplots(1,1, figsize=(7,5))
      ax = xAx[i]
      counter = 0
      #for solutionfile in os.listdir(workDir):
      #  if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
      #    assert counter==0,"Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

      #    hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),VARKEY)
      #    print "Interpolating"
      #    vq = ci.interpolate(hr,vtx,wts)

      #    print "axi mean"
      #    xQtemp = []
      #    for it in range(rest):
      #      imin = it * resx*resr
      #      imax = (it + 1) * resx*resr
      #      xQtemp.append(vq[imin:imax])
      #    xQtemp = np.array(xQtemp)
      #    mean = np.zeros(resx*resr)


      #    for ip in range(resx*resr):
      #      mean[ip] = np.sum(xQtemp[:,ip])/rest

      #    Q = np.reshape(mean,(resr,resx))

      #    Q[np.isnan(Q)] = 0

      profile = []
      for j in range(resr):
        interp = interp1d(xq, Q[j,:], kind='cubic')
        profile.append(interp(normAxialDistance*D))
      ax.plot(rq, profile, label=label)
      ax.legend(fontsize=12)
      #counter +=1
      #  else:
      #    continue




def MasterPlotProfile(workDir, D, variable, label, xAx, xNX, interpRelPath, solRelPath=None):
  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'av_sol'

  # Save images to separate folder
  IMG_OUT = workDir
  
  assert len(xAx)==len(xNX)+1 , 'The number of inputed figs must be the same as the number of axial positions'
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



  ax = xAx[0]
  counter = 0
  for solutionfile in os.listdir(workDir):
    if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
      assert counter == 0, "Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

      hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir), solutionfile), VARKEY)
      print "Interpolating"
      vq = ci.interpolate(hr, vtx, wts)

      print "axi mean"
      xQtemp = []
      for it in range(rest):
        imin = it * resx * resr
        imax = (it + 1) * resx * resr
        xQtemp.append(vq[imin:imax])
      xQtemp = np.array(xQtemp)
      mean = np.zeros(resx * resr)

      for ip in range(resx * resr):
        mean[ip] = np.sum(xQtemp[:, ip]) / rest

      Q = np.reshape(mean, (resr, resx))

      Q[np.isnan(Q)] = 0

      ax.plot(xq, Q[0,:], label=label)
      ax.legend(fontsize=12)
      counter += 1
    else:
      continue

  for i in range(1,len(xNX)+1):
    normAxialDistance = xNX[i-1]
    #fig, ax = plt.subplots(1,1, figsize=(7,5))
    ax = xAx[i]
    counter = 0
    #for solutionfile in os.listdir(workDir):
    #  if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX):
    #    assert counter==0,"Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

    #    hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir),solutionfile),VARKEY)
    #    print "Interpolating"
    #    vq = ci.interpolate(hr,vtx,wts)

    #    print "axi mean"
    #    xQtemp = []
    #    for it in range(rest):
    #      imin = it * resx*resr
    #      imax = (it + 1) * resx*resr
    #      xQtemp.append(vq[imin:imax])
    #    xQtemp = np.array(xQtemp)
    #    mean = np.zeros(resx*resr)


    #    for ip in range(resx*resr):
    #      mean[ip] = np.sum(xQtemp[:,ip])/rest

    #    Q = np.reshape(mean,(resr,resx))

    #    Q[np.isnan(Q)] = 0

    profile = []
    for j in range(resr):
      interp = interp1d(xq, Q[j,:], kind='cubic')
      profile.append(interp(normAxialDistance*D))

    ax.plot(rq, profile, label=label)
    ax.legend(fontsize=12)
    #    counter +=1
    #  else:
    #    continue

  
  



def PlotNormalizedAxialVelocityProfileOnAxis(workDir, D, label, xAx, interpRelPath, isFirst, solRelPath=None, virtualSource=0.0, uCoflow=3.5):
  # Set general path to solution files and interpolation base file
  FILEPREFIX = 'av_sol'

  # Save images to separate folder
  IMG_OUT = workDir
 
  xNX=[0,0,0] 
  assert len(xAx)==len(xNX), 'The number of inputed figs must be the same as the number of axial positions'
  # Keyname of variable to integrate (as specified in solution h5-file)
  VARKEY = 'Average/u'


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
      assert counter == 0, "Too much solution files in this directory (" + workDir + "). PANIIIIIC !"

      hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(workDir), solutionfile), VARKEY)
      print "Interpolating"
      vq = ci.interpolate(hr, vtx, wts)

      print "axi mean"
      xQtemp = []
      for it in range(rest):
        imin = it * resx * resr
        imax = (it + 1) * resx * resr
        xQtemp.append(vq[imin:imax])
      xQtemp = np.array(xQtemp)
      mean = np.zeros(resx * resr)

      for ip in range(resx * resr):
        mean[ip] = np.sum(xQtemp[:, ip]) / rest

      Q = np.reshape(mean, (resr, resx))

      Q[np.isnan(Q)] = 0

      #Substracting coflow velocity
      Q = Q - uCoflow

      U0x = Q[0,:]
      U0  = U0x[0]
      
      profile = U0 / U0x

      x0 = virtualSource
      abscissa = xq - x0

      # Profils axiaux
      xAx[0].plot(abscissa, profile, label=label)

      if isFirst:
        B   = 5.8 # cf Pope p.101 For round jets of uniform density
        xx  = np.linspace(abscissa.min(), abscissa.max(), 100)
        yy  = xx / (B * D)
        xAx[0].plot(xx, yy, color='k', label='Theory - B = 5.8')

      xAx[0].set_xlabel('$x - x_0$')
      xAx[0].set_ylabel('$U_0/U_0(x)$')

      xAx[0].legend(fontsize=12)




      # Elargissement du jet
      r12 = np.zeros_like(U0x)

      rProfile = []
      for i in range(resx):
        interp = interp1d(Q[:,i], rq, kind='cubic')
        rProfile.append( interp( U0x[i]/2.0 ) )


      xAx[1].plot(abscissa, rProfile, label=label)

      if isFirst:
        S   = 0.094 # cf Pope p.101 For round jets of uniform density
        xx  = np.linspace(abscissa.min(), abscissa.max(), 100)
        yy  = xx * S
        xAx[1].plot(xx, yy, color='k', label='Theory - S = 0.094')
        xAx[1].set_ylim([yy.min(), yy.max()])

      xAx[1].set_xlabel('$x - x_0$')
      xAx[1].set_ylabel('$r_{1/2}(x)$')

      xAx[1].legend(fontsize=12)



      # Profils autosimilaires

      xXAdim = np.array([8,10,12,14, 16, 18, 20, 22])
      xX = xXAdim * D

      for i in range(len(xX)):
        axialDistance = xX[i]
        # Determination du profil a tracer
        profile = []
        for j in range(resr):
          interp = interp1d(xq, Q[j,:], kind='cubic')
          profile.append(interp(axialDistance))
        profile = np.array(profile)
        # Determination de la vmax du profil
        U0x = profile[0]
        # Determination du demi-rayon 
        interp2 = interp1d(profile, rq, kind='cubic')
        r12 = interp2( U0x/2.0  )

        radim = rq / r12
        profadim = profile / U0x

        xAx[2].plot(radim, profadim, label=label + ' - r/d = ' + str(xXAdim[i]))

      xAx[2].legend(fontsize=12)
      xAx[2].set_xlim([0.0,3.0])
      xAx[2].set_xlabel('$r / r_{1/2}(x)$')
      xAx[2].set_ylabel('$U / U_0(x)$')

      xAx[0].vlines(xX, 0.0, 7)
      xAx[1].vlines(xX, 0.002, 0.015)

      counter += 1
    else:
      continue
