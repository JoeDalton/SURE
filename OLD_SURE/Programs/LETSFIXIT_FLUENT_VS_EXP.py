#!/usr/bin/env python2.7
import argparse
import Link
import matplotlib.pyplot as plt
from AVBP_Tools.AxiFrom3D.PlotProfile import PlotProfile
import numpy as np
#import Utils.PrettyPlots as pp
from Utils.I_O import ErrorMsg, LoadCSV 
import os
from scipy.interpolate import interp1d
from scipy.stats import linregress


# Parse arguments
parser = argparse.ArgumentParser(description="Explicit file name")
parser.add_argument('--fd'   , required=True, nargs='+', type=str , help = "List of dirs for fluent")
parser.add_argument('--flab' , required=True, nargs='+', type=str , help = "List of labels corresponding to the fluent dirs")
parser.add_argument('--ed'   , required=True, type=str , help = "Exp dir")

args  = parser.parse_args()
xfd   = args.fd
xflab = args.flab
ed    = args.ed

D = 0.00457


# Read data
lXExp     = ["01","05","10","15","25"]
lExpArray = []
for pos in lXExp:
  fileName  = ed + "/T_" + pos + ".csv"
  x,y       = LoadCSV(fileName)
  lExpArray.append([x,y])

lXFl      = ["01", "05", "10", "15", "20", "25", "30", "35", "40"]
llFluentArray = []
for wd in xfd:
  wdArray = []
  for pos in lXFl:
    fileName  = wd + "/T_" + pos + ".dat"
    if os.path.isfile(fileName):
      x,y       = LoadCSV(fileName, delimiter=" ", nHeaderLines=4)
      x         = x / D
      sortorder = np.argsort(x)
      x = x[sortorder]
      y = y[sortorder]
      wdArray.append([x,y])
    else:
      wdArray.append([])
  llFluentArray.append(wdArray)




# Extract axial profile
xExpAxial = []
yExpAxial = []
for i in range(len(lExpArray)):
  p = lExpArray[i]
  xExpAxial.append(float(lXExp[i]))
  yExpAxial.append(p[1][0])


xXFlAxial = []
xYFlAxial = []
for array in llFluentArray:
  xFlAxial = []
  yFlAxial = []
  for i in range(len(array)):
    p = array[i]
    xFlAxial.append(float(lXFl[i]))
    yFlAxial.append(p[1][0])
  xXFlAxial.append(xFlAxial)
  xYFlAxial.append(yFlAxial)


# Plot axial profile
fig = plt.figure(figsize=(8,7))
ax = plt.subplot(111)

ax.set_aspect('auto')
ax.set_ylabel("Axial Temperature [K]")
ax.set_xlabel("$x/d$ [-]")

ax.scatter(xExpAxial,yExpAxial, label="Exp", edgecolors='k',facecolors='none', marker='o')
for i in range(len(xflab)):
  ax.plot(xXFlAxial[i], xYFlAxial[i], label=xflab[i])

#Linear regression over Exp data and shift for num data
regr = linregress(xExpAxial, yExpAxial)
slope, intercept = regr.slope, regr.intercept
newX = np.linspace(np.min(xExpAxial), np.max(xExpAxial), 1000)
ax.plot(newX, slope*newX + intercept, label='Regression for Exp')
ax.plot(newX+20, slope*(newX) + intercept, label='Shift for k-w-std')



ax.legend()

plt.title("Axial temperature")

plt.savefig('Axial.pdf', bbox_inches="tight", pad_inches = 0)
plt.cla()
plt.clf()
plt.close()


# Extract jet width (w.r.t. Temperature)
maxTemp = 1190.0


xExpWidth = []
yExpWidth = []
for i in range(len(lExpArray)):
  centerTemp  = yExpAxial[i]
  halfTemp    = (maxTemp + centerTemp) / 2.0
  xp          = lExpArray[i][0]
  yp          = lExpArray[i][1]
  maxI = np.argmax(yp)
  yp[maxI-1] = yp[maxI-1] + 0.01 #DEBUG...
  try:
    fwidth      = interp1d(yp[:maxI],xp[:maxI], kind='linear')
  except:
    print yp[:maxI]
    quit()
  try:
    width       = fwidth(halfTemp)
  except:
    print halfTemp
    print yp[:maxI]
    quit()
  xExpWidth.append(float(lXExp[i]))
  yExpWidth.append(width)









xXFlWidth = []
xYFlWidth = []
for j in range(len(llFluentArray)):
  array = llFluentArray[j]
  xFlWidth = []
  yFlWidth = []
  for i in range(len(array)):
    centerTemp  = xYFlAxial[j][i]
    xp = array[i][0]
    yp = array[i][1]

    maxI = np.argmax(yp)
    yp[maxI-1] = yp[maxI-1] + 0.01 #DEBUG...
    try:
      fwidth      = interp1d(yp[:maxI],xp[:maxI], kind='linear')
    except:
      print yp[:maxI]
      quit()
    try:
      width       = fwidth(halfTemp)
    except:
      print halfTemp
      print yp[:maxI]
      quit()
    xFlWidth.append(float(lXFl[i]))
    yFlWidth.append(width)
    
  xXFlWidth.append(xFlWidth)
  xYFlWidth.append(yFlWidth)


# Plot axial profile
fig = plt.figure(figsize=(8,7))
ax = plt.subplot(111)

ax.set_aspect('auto')
ax.set_ylabel("$r_{1/2}/d$ [-]")
ax.set_xlabel("$x/d$ [-]")

ax.scatter(xExpWidth,yExpWidth, label="Exp", edgecolors='k',facecolors='none', marker='o')
for i in range(len(xflab)):
  ax.plot(xXFlWidth[i], xYFlWidth[i], label=xflab[i])

#Linear regression over Exp data and shift for num data
regr = linregress(xExpWidth, yExpWidth)
slope, intercept = regr.slope, regr.intercept
newX = np.linspace(0, 40, 1000)
ax.plot(newX, slope*newX + intercept, label='Regression for Exp')
ax.plot(newX, slope*(newX) + intercept+0.4, label='Shift for k-w-std')

ax.legend()

plt.title("Jet half-width w.r.t temperature")

#plt.show()
plt.savefig('Width.pdf', bbox_inches="tight", pad_inches = 0)
fig = None
