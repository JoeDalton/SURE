import os
import numpy as np
import InterpolateFromAVBP as ci
import common_IO as cio
import math






def GenInterpFile_Slice(workDir, meshRelPath, interpRelPath, xlim, ylim):
  # Dimension constant. DO NOT CHANGE! (unless you know what you're doing)
  d=2

  # Set resolution of interpolation mesh
  resx = 315
  resy = 500  # Double what I take for the axi mesh because I interpolate both sides

  # Read data
  print(workDir)
  print(meshRelPath)
  x, y, z = cio.read_avbp_mesh_coordinates(os.path.join(workDir + '/' + meshRelPath))

  print "Stacking coordinates"
  #xyz = np.vstack((x,y,z)).T
  xy = np.vstack((x,y)).T

  # Generate interpolation mesh
  print "Generate interpolation mesh"
  xq = np.linspace(xlim[0],xlim[1],resx)
  yq = np.linspace(ylim[0],ylim[1],resy)
  #zq = np.linspace(zlim[0],zlim[1],resz)

  print "Initial mesh contains %i nodes" % (len(x))
  print "Interpolation mesh contains %i nodes" % (resx*resy)

  # Create coordinates of interpolation mesh
  # Note: Unlike previous implementation (using meshgrid),
  #       x,y,z vectors are created, no matrices!
  uv = np.vstack(np.meshgrid(xq, yq)).reshape(2,-1).T

  # Compute the weights for each vertex of the interpolation mesh
  print "Triangulating and calculating weights. This may take a while..."
  vtx, wts = ci.interp_weights(xy, uv, d)

  # Dump vertices and weights
  header = np.array([['resx', 'resy', 'xlim', 'ylim', 'isAxiSym'], [resx, resy, xlim, ylim, True]], dtype=object)
  print(workDir+'/'+interpRelPath)
  print(len(vtx))
  print(len(wts))
  print(len(header))
  cio.save_interp_base(workDir+'/'+interpRelPath, vtx=vtx, wts=wts, header=header, has_voln=False, voln=None)

  print "Done!"
