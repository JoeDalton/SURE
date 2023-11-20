import os
import numpy as np
import InterpolateFromAVBP as ci
import common_IO as cio
import math






def GenInterpFile_Cart(workDir, meshRelPath, interpRelPath, xlim, ylim, zlim):
  # Dimension constant. DO NOT CHANGE! (unless you know what you're doing)
  d=3

  # Set resolution of interpolation mesh
  #resx = 300
  #resr = 250
  #rest = 20
  resx = 315
  resy = 500  # Double what I take for the axi mesh because I interpolate both sides
  resz = 3    # Only 3 : this way I can take the 2nd row to get the slice at z = 0

  # Read data
  print(workDir)
  print(meshRelPath)
  x, y, z = cio.read_avbp_mesh_coordinates(os.path.join(workDir + '/' + meshRelPath))

  print "Stacking coordinates"
  xyz = np.vstack((x,y,z)).T

  # Generate interpolation mesh
  print "Generate interpolation mesh"
  xq = np.linspace(xlim[0],xlim[1],resx)
  yq = np.linspace(ylim[0],ylim[1],resy)
  zq = np.linspace(zlim[0],zlim[1],resz)

  print "Initial mesh contains %i nodes" % (len(x))
  print "Interpolation mesh contains %i nodes" % (resx*resy*resz)

  # Create coordinates of interpolation mesh
  # Note: Unlike previous implementation (using meshgrid),
  #       x,y,z vectors are created, no matrices!
  uvw = np.vstack(np.meshgrid(xq, yq, zq)).reshape(3,-1).T

  # Compute the weights for each vertex of the interpolation mesh
  print "Triangulating and calculating weights. This may take a while..."
  vtx, wts = ci.interp_weights(xyz, uvw, d)

  # Dump vertices and weights
  header = np.array([['resx', 'resy', 'resz', 'xlim', 'ylim', 'zlim', 'isAxiSym'], [resx, resy, resz, xlim, ylim, zlim, True]], dtype=object)
  print(workDir+'/'+interpRelPath)
  print(len(vtx))
  print(len(wts))
  print(len(header))
  cio.save_interp_base(workDir+'/'+interpRelPath, vtx=vtx, wts=wts, header=header, has_voln=False, voln=None)

  print "Done!"
