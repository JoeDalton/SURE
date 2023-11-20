#!/usr/bin/env python2.7

import Link
import numpy                as np
import openturns            as ot
import matplotlib.pyplot    as plt
from Sampling.MC            import *
from Sampling.QMC           import *
from Sampling.RQMC          import *
from Sampling.Deterministic import *


################## Cubature
dim                 = 3
#rule                = 'Clenshaw-Curtis'
rule                = 'SecondFejer'
distributionType    = 'Uniform' #'LogNormal'
xDistributionParam  = [0.0,1.0]
level               = 3 #4
xRule               = [rule] * dim
xDistributionType   = [distributionType] * dim
xxDistributionParam = [xDistributionParam] * dim
xLevel              = [level] * dim
#xLevel              = [level, level]

# Full
#CubDrawer = CubatureSampling(None, xRule, xDistributionType, xxDistributionParam, xLevel=xLevel, isSparse = False, verbosity = True)
#CubDrawer = CubatureSampling(None, xRule, xDistributionType, xxDistributionParam, maxLevel=level, isSparse = False, verbosity = True)

#Smolyak
CubDrawer = CubatureSampling(None, xRule, xDistributionType, xxDistributionParam, xLevel=xLevel, isSparse = True, verbosity = -1)
#CubDrawer = CubatureSampling(None, xRule, xDistributionType, xxDistributionParam, maxLevel=level, isSparse = True, verbosity = True)
CubDrawer.Draw()

x = np.zeros(np.size(CubDrawer.xPoint, 0))
x[:] = CubDrawer.xPoint[:,0]
y = np.zeros(np.size(CubDrawer.xPoint, 0))
y[:] = CubDrawer.xPoint[:,1]
z = np.zeros(np.size(CubDrawer.xPoint, 0))
z[:] = CubDrawer.xPoint[:,2]


#--------------------------------------------------------------------- Plot
# Cubature
#fig = plt.figure(figsize=(14,10))
#plt.scatter(x, y, label='Cubature')
#plt.legend(loc='best')
#plt.show()


from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
# Create a figure and a 3D Axes
fig = plt.figure()
ax = Axes3D(fig)

# Create an init function and the animate functions.
# Both are explained in the tutorial. Since we are changing
# the the elevation and azimuth and no objects are really
# changed on the plot we don't have to return anything from
# the init and animate function. (return value is explained
# in the tutorial.
def init():
    ax.scatter(x, y, z, marker='o', s=20, c="k", alpha=0.6)
    return fig,

#def animate(i):
#    ax.view_init(elev=10., azim=i)
#    return fig,
#
#
#
## Animate
#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=360, interval=20, blit=True)
init()
plt.show()

# Save
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
