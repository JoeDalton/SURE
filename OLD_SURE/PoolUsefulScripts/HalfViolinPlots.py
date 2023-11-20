
import numpy as np
import os
import subprocess
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter as cc
import matplotlib as mpl
font = {'family' : 'sans',
    'serif': 'helvet',
    'weight' : 'bold',
    'size'   : 20}

mpl.rc('font', **font)
mpl.rc('text', usetex=True)

precursorList = [1015,1030,1045,1060]
tList = precursorList


resultFileName = "chemOnly.h5"

data2019  =[]
data2019N =[]
for T in tList:
  tFolderName	= 'Konnov_2019/T' + str(T) + '/results'
  tFile	= tFolderName + "/" + resultFileName
  res		= h5py.File(tFile, 'r')
  Raw		= res["Data/QoI"][:]
  Raw		= np.log10(Raw[Raw!=0.0])
  res.close()
  tFolderName	= 'Konnov_2019_naked_Unc/T' + str(T) + '/results'
  tFile	= tFolderName + "/" + resultFileName
  res		= h5py.File(tFile, 'r')
  nak		= np.log10(res["Data/QoI"][:])
  res.close()
  data2019.append(Raw)
  data2019N.append(nak)


data1 = np.array(data2019).T #(np.random.normal(0, 1, size=10000), np.random.normal(0, 2, size=10000))
data2 = np.array(data2019N).T #(np.random.normal(1, 1, size=10000), np.random.normal(1, 2, size=10000))

fig, ax = plt.subplots(figsize=(7, 5))

v1 = ax.violinplot(data1, points=1000, positions=precursorList,
               showmeans=False, showextrema=False, showmedians=False, widths=5)
for b in v1['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further right than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
    b.set_facecolor('#D43F3A')
    b.set_edgecolor('black')
    b.set_alpha(1)

v2 = ax.violinplot(data2, points=1000, positions=precursorList, 
               showmeans=False, showextrema=False, showmedians=False, widths=5)

for b in v2['bodies']:
    # get the center
    m = np.mean(b.get_paths()[0].vertices[:, 0])
    # modify the paths to not go further left than the center
    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
    b.set_color('b')
    b.set_facecolor('#4169E1')
    b.set_edgecolor('black')
    b.set_alpha(1)

ax.legend([v1['bodies'][0],v2['bodies'][0]],['Vanilla Konnov (2019)', 'Arranged Konnov (2019)'])

ax.set_xlabel(r'$T_{co-flow}$ [K]')
ax.set_ylabel(r'log$_{10}$($\tau$) [-]')
#ax.grid(True)
for item in [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels():
  item.set_fontsize(20)

plt.savefig('K2019vsK2019N.pdf', bbox_inches='tight')
