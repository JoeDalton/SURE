import os
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
#import Utils.PrettyPlots as pp
import AVBP_Tools.AxiFrom3D.InterpolateFromAVBP as ci
import AVBP_Tools.AxiFrom3D.common_IO as cio




def Scatterplot(sol, mesh, title, xV):
  
  if 'x' in xV or 'y' in xV or 'z' in xV:
    # Read mesh
    x, y, z = cio.read_avbp_mesh_coordinates(mesh)
  
  if xV[0] == 'x':
    var1 = x
  elif xV[0] == 'y':
    var1 = y
  elif xV[0] == 'z':
    var1 = z
  else:
    var1 = cio.read_avbp_solution(sol,xV[0])

  if len(xV) == 1:
    #plot histogram
    plt.hist(var1, density=True, color='k')
    plt.show()
    quit()

  ####################################################################
  if xV[1] == 'x':
    var2 = x
  elif xV[1] == 'y':
    var2 = y
  elif xV[1] == 'z':
    var2 = z
  else:
    var2 = cio.read_avbp_solution(sol,xV[1])

  if len(xV) == 2:
    fig = plt.figure(1, figsize=(15,10))
    ax0 = fig.add_subplot(222)
    ax0.scatter(var1, var2, marker=',', s=1, c='k')
    ax0.set_xlabel(xV[0])#, fontsize=16)
    ax0.set_ylabel(xV[1])#, fontsize=16)
    ax1 = fig.add_subplot(221)
    ax1.hist(var2, density=True, color='k', orientation='horizontal')
    ax1.set_xlabel('pdf')#, fontsize=16)
    ax1.set_ylabel(xV[1])#, fontsize=16)
    ax1.invert_xaxis()
    ax2 = fig.add_subplot(224)
    ax2.hist(var1, density=True, color='k')
    ax2.invert_yaxis()
    ax2.set_xlabel(xV[0])#, fontsize=16)
    ax2.set_ylabel('pdf')#, fontsize=16)
    plt.show()
    quit()

  ####################################################################
  if xV[2] == 'x':
    var3 = x
  elif xV[2] == 'y':
    var3 = y
  elif xV[2] == 'z':
    var3 = z
  else:
    var3 = cio.read_avbp_solution(sol,xV[2])

  #fig, ax = plt.subplots(1, 1, figsize=(15, 10))

  cmap='inferno'


  #ax.scatter(var1, var2, marker=',', s=1, c=var3, cmap=cmap)
  #norm = mpl.colors.Normalize(vmin=var3.min(),vmax=var3.max())
  #sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
  #sm.set_array([])
  #cbar = plt.colorbar(sm)
  #cbar.set_label(xV[2])
  #plt.title(title)
  #plt.xlabel(xV[0])#, fontsize=16)
  #plt.ylabel(xV[1])#, fontsize=16)



  fig = plt.figure(1, figsize=(15,10))
  ax0 = fig.add_subplot(222)
  ax0.scatter(var1, var2, marker=',', s=1, c=var3, cmap=cmap)
  #ax0.scatter(var1, var2, marker=',', s=1, c='k')
  ax0.set_xlabel(xV[0])#, fontsize=16)
  ax0.set_ylabel(xV[1])#, fontsize=16)
  norm = mpl.colors.Normalize(vmin=var3.min(),vmax=var3.max())
  sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
  sm.set_array([])
  cbar = plt.colorbar(sm)
  cbar.set_label(xV[2])

  ax1 = fig.add_subplot(221)
  ax1.hist(var2, density=True, color='k', orientation='horizontal')
  ax1.set_xlabel('pdf')#, fontsize=16)
  ax1.set_ylabel(xV[1])#, fontsize=16)
  ax1.invert_xaxis()
  ax2 = fig.add_subplot(224)
  ax2.hist(var1, density=True, color='k')
  ax2.invert_yaxis()
  ax2.set_xlabel(xV[0])#, fontsize=16)
  ax2.set_ylabel('pdf')#, fontsize=16)


  ax3 = fig.add_subplot(223)
  ax3.hist(var3, density=True, color='grey')
  ax3.set_xlabel(xV[2])#, fontsize=16)
  ax3.set_ylabel('pdf')#, fontsize=16)



  plt.show()
  quit()

  #plt.savefig(title + '.png', bbox_inches='tight', pad_inches=0)
  


