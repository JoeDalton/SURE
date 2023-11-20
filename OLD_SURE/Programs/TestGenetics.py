#!/usr/bin/env python2.7

import Link
import numpy                as np
import matplotlib.pyplot    as plt
import random
from scipy.interpolate                          import RectBivariateSpline
from Numerics.Optimisation.Genetics.Genetics    import *
from Numerics.Optimisation.Genetics.Population  import *
from Numerics.Optimisation.Genetics.Eval        import MultiEval

# ===============================================
# =========    Function definition    ===========
# ===============================================
xMin  = [0.0,0.0]
xMax  = [1.0,1.0]
#######################
##  "Random" Spline  ##
#######################
resol = 10
random.seed(0)
x0 = np.linspace(xMin[0],xMax[0],resol)
y0 = np.linspace(xMin[1],xMax[1],resol)
z0 = np.zeros((resol,resol))
for i in range(resol):
  for j in range(resol):
    z0[i,j] = random.random()
MySpl = RectBivariateSpline(x0, y0, z0, s=0)
#######################
##    Tilted plane   ##
#######################
def MyPlane(x,y):
  z = x + y
  return z

MyFunc = MySpl
#MyFunc = MyPlane

# ===============================================
# =========     Optimisation type     ===========
# ===============================================
def Opti1(vect):
  # vect must have 2 components
  result = MyFunc(vect[0], vect[1])
  xSign = [-1]
  xCrit = []
  xGrade = [result]
  return MultiEval(xGrade, xCrit, xSign=xSign)

threshold = 0.7
tolerance = 1.0e-3
def Opti2(vect):
  # vect must have 2 components
  # grade 1: simulation must have a certain threshold value
  grade1 = abs(MyFunc(vect[0], vect[1]) - threshold)
  # grade 2: x must be as small as possible
  grade2 = vect[0]
  xSign = [1, 1]
  xCrit = [tolerance]
  xGrade = [grade1, grade2]
  return MultiEval(xGrade, xCrit, xSign=xSign)

EvalFunction = Opti1
# ===============================================
# ==============    Parameters    ===============
# ===============================================
outLevel        = 2
#preSet = "Hold The Line" # "I'm Too Young To Die" | "Hurt Me Plenty" | "Nightmare" | "Hold The Line"
preSet = "Hurt Me Plenty"

# ===============================================
# ==============  Initialization  ===============
# ===============================================
myPopulation  = Population(xMin, xMax, EvalFunction, preSet=preSet)
myAlgorithm   = GeneticAlgorithm(myPopulation, outLevel=outLevel, preSet=preSet)


# ===============================================
# ==============       Run        ===============
# ===============================================
#myAlgorithm.Run(parallel=True)
myAlgorithm.Run(parallel=False)


# ===============================================
# ==============     Results      ===============
# ===============================================
bestXY    = myPopulation.xRat[0].xGen
bestValue = EvalFunction(bestXY).xGrade[0]



# ===============================================
# ==============       Plot       ===============
# ===============================================




from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

nx  = 100
ny  = 100
x   = np.linspace(xMin[0], xMax[0], resol*10)
y   = np.linspace(xMin[1], xMax[1], resol*10)
X,Y = np.meshgrid(x,y)
Z   = np.zeros_like(X)
for i in range(nx):
  for j in range(ny):
    Z[i,j] = MyFunc(X[i,j],Y[i,j])

norm = plt.Normalize(Z.min(), Z.max())
colors=cm.hot(norm(Z))
rcount, ccount, _ = colors.shape

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(111,projection='3d')
ax.plot_wireframe(X, Y, Z, cmap=cm.coolwarm, alpha = 0.3)
ax.scatter([bestXY[0]], [bestXY[1]], [MyFunc(bestXY[0],bestXY[1])], marker = 'o', color='k')

#cset = ax.contour(X, Y, Z, 0, zdir='z', offset=0.5, cmap=cm.coolwarm)
cset = ax.contour(X, Y, Z, threshold, zdir='z', cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, 20, zdir='x', offset=-4, cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, 20, zdir='y', offset=1.1, cmap=cm.coolwarm)

ax.set_xlabel(r'x', labelpad=15)
ax.set_ylabel(r'y', labelpad=15)
ax.set_zlabel(r'z', labelpad=15)
plt.show()
#plt.savefig(ModelFile, bbox_inches='tight')
#plt.cla()
#plt.clf()
#plt.close()
