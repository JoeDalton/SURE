#!/usr/bin/env python2.7

import Link
import numpy                as np
import matplotlib.pyplot    as plt
import openturns            as ot
from Sampling.MC            import *
from Sampling.QMC           import *
from Sampling.RQMC          import *
from Sampling.Deterministic import *


#------------------------------------------------------------ Preparation of samples
nSample      = 50
firstSample  = 0
nQMCDrawer   = 10
distribution = ot.ComposedDistribution([ot.Uniform(0.0, 1.0)] * 2)
verb = -1


############### QMC
QMCDrawer  = QMC(None, distribution, nSample, firstSample=firstSample, sequence='Sobol', verbosity=verb)
s2 = QMCDrawer.Draw()
y2 = np.zeros(nSample)
x2 = np.zeros(nSample)
x2 = s2[:, 0]
y2 = s2[:, 1]

# QMC
fig = plt.figure(figsize=(7,5))
plt.scatter(x2, y2, facecolors='none', edgecolors='k')
plt.savefig('sobol.pdf', bbox_inches="tight")

