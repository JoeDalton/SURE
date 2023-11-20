import h5py as h5
import os
import numpy as np
import matplotlib
#matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc
import custom_interpolation as ci
import common_IO as cio
import math
import csv


##########################################################################################

# 	USER INPUT

##########################################################################################

# Set general path to solution files and interpolation base file
pathname = '/Users/lavabreg/Desktop/Thesis/TempStudy_K2019N/'
FNAME_INTERP_BASE = '/Users/lavabreg/Desktop/Thesis/TempStudy_K2019N/toto.npz'
FILEPREFIX = 'av_sol'

# Keyname of variable to integrate (as specified in solution h5-file)
VARKEY = 'Average/c_out'

# Integration axis (default = 2 corresponding to z-direction)
INT_AXIS = 2

# Activate voln correction, if desired AND available
activate_voln = False

# Nozzle diameter
D=0.00457

# Axis subdivisions
YSSIZE=1
XSSIZE=1


xstep = 10
rstep = 10

# Subdivisions of colormap
CMAP_DIV = 50

# Colormap legend
CMAP_LEGEND = 'Normalized Progress Variable [-]'

# Save with fancy colormap as PDF
MAKE_PDF = False

# Save images to separate folder
IMG_OUT = pathname

# Prepend identifier for images (underscore is automatically inserted)
IM_FNAME = 'axi_mean'














##########################################################################################



# Read data
interp_base = cio.load_interp_base(FNAME_INTERP_BASE)
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
tq = np.linspace(tlim[0],tlim[1],rest+1)[:-1]
#tq = np.linspace(0.0,360.0,rest+1)[:-1]

print "Interpolation mesh contains %i cells" % (resx*resr*rest)

# Loop through all solution files and apply interpolation on coarse grid

qoiList=[]

for solutionfile in os.listdir(pathname):
  print(solutionfile)
  if solutionfile.endswith(".h5") and solutionfile.startswith(FILEPREFIX): 
    hr = cio.read_avbp_solution(os.path.join(os.path.expanduser(pathname),solutionfile),VARKEY)
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
    #        rc = {"figure.subplot.left"    : 0.1, 
    #              "figure.subplot.right"   : 0.9,
    #              "figure.subplot.bottom"  : 0.1,
    #              "figure.subplot.top"     : 0.9 }
    if MAKE_PDF:
      pass
    else:	
      fig, ax = plt.subplots(1,1, figsize=(7,7))
      #font = {'family' : 'sans',
      #    'serif': 'helvet',
      #    'weight' : 'bold',
      #    'size'   : 16}
      #
      #plt.rc('font', **font)
      #plt.rc('text', usetex=True)
      cmap = 'hot'
      norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
      ax = plt.gca()
      img = ax.imshow(np.rot90(Q,1), interpolation='bilinear', cmap=cmap, norm=norm, extent= [rq.min(), rq.max(), xq.min(), xq.max()])
      plt.gca().set_aspect('equal', adjustable='box')

      #plt.contourf(rq, xq, np.rot90(Q,1), [0.1])
      myContour = plt.contour(rq, xq, np.flipud(np.rot90(Q,1)), [0.1], colors='w', linestyles='--', linewidths=1.5)

      cbar = plt.colorbar(img)
      cbar.set_label(CMAP_LEGEND,size=16)
      cbar.ax.tick_params(labelsize=12)
      
      xticks = np.arange(rlim[0], rlim[1], rstep*D)
      xtickslabels = np.round(xticks/D)
      plt.xticks(xticks, xtickslabels, fontsize=12)
      
      yticks = np.arange(xlim[0], xlim[1], xstep*D)
      ytickslabels = np.round(yticks/D)
      plt.yticks(yticks, ytickslabels, fontsize=12)
      
      plt.xlabel('r/D',fontsize=16)
      plt.ylabel('x/D',fontsize=16) 


      """
      Determine lift-off height
      """
      qoi=min(myContour.allsegs[0][0][:,1])/D
      print(qoi)
     
      plt.title('H/D = ' + str(qoi)) 
      plt.hlines(y=qoi*D, xmin=0, xmax=0.1, color='w', linestyle=":")

      fname=os.path.split(solutionfile)[1].split('.')[0]
      plt.savefig(IMG_OUT + '/' + IM_FNAME + '_%s.pdf' % fname, bbox_inches="tight",pad_inches = 0)

      qoiList.append(qoi)

  else:
    continue


qoiList.sort()
qoiList.reverse()
#TList = [1015,1030,1045,1060]
TList = [1015, 1030, 1045]
rows = np.array([TList, qoiList]).T
header=['T', 'QOI']

csvfile=pathname+'K2019Results.csv'
with open(csvfile, 'w') as csvf:
  csvw = csv.writer(csvf)
  csvw.writerow(header)
  csvw.writerow(header)
  csvw.writerow(header)
  csvw.writerow(header)
  csvw.writerow(header)
  csvw.writerow(header)
  csvw.writerows(rows)
##fig, ax = plt.subplots(1,1, figsize=(7,5))
##plt.plot(TList, qoiList, color='k', marker='p')
##plt.show()

print "Done!" 
