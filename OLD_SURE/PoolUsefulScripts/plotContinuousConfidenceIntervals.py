#!/usr/bin/env python2.7

import numpy as np
import os
import subprocess
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter as cc

def plot_mean_and_CI(mean, lb, ub, tList, color_mean=None, color_shading=None):
    abscissa = np.zeros_like(tList)
    abscissa[:] = 1000 / tList[:]
    # plot the shaded range of the confidence intervals
    plt.fill_between(abscissa, ub, lb,
                     color=color_shading, alpha=.5)
    # plot the mean on top
    plt.plot(abscissa, mean, color_mean, marker = '+')
    plt.yscale('log')

# ----- Preventing Macs from going to sleep too early... -----
if os.environ["COMMCOMB_PLATFORM"]=="MAC_LIGHT" :
    coffee = subprocess.Popen('caffeinate')

precursorList = np.linspace(0.95, 1.2, 30)
tList = 1000/precursorList

tList[:] = np.around(tList[:], 1)

data = []

# ----- Reading solution files -----
for T in tList:
    mainFolderName = 'T_' + str(T)
    print mainFolderName
    counter = 0
    folderName = mainFolderName + '/sample_0'
    IDT_List=[]
    while os.path.isdir(folderName):
	filePath = './' + folderName + '/Solution.h5'
	try:
	    sol = h5py.File(filePath, 'r')
	    IDT = float(sol['Header/Main_IgnitionDelay'][0])
	    IDT_List.append(IDT)
	    sol.close()
	except:
	    print('Invalid solution file.')
	counter += 1
	folderName = mainFolderName + '/sample_' + str(counter)
    
    if counter==0:
	print "no matching folder : " + mainFolderName

    data.append(IDT_List)
# ----- Compute mean and confidence intervals -----
mean = np.zeros([len(data)])
lb = np.zeros([len(data)])
ub = np.zeros([len(data)])

for i in range(len(data)):
    mean[i] = np.average(data[i])
    lb[i] = np.percentile(data[i], 0.5)
    ub[i] = np.percentile(data[i], 99.5)

fig = plt.figure(1, figsize=(7, 5))
plot_mean_and_CI(mean, ub, lb, tList, color_mean='k', color_shading='k')
#plt.boxplot(IDT, positions = 1000/tList, widths=0.01, whis=[0.5,99.5])
plt.title('Auto-ignition delay as a function of initial temperature')
plt.xlabel('1000/T')
plt.ylabel('IDT (s)')
plt.show()

# ----- Termination -----
if os.environ["COMMCOMB_PLATFORM"]=="MAC_LIGHT":
    coffee.kill()
