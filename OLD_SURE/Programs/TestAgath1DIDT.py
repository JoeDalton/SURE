#!/usr/bin/env python2.7

import Link
from Interface.Results  import AGATH_1D_IDT
from Utils.I_O          import prt

solutDir    = "."
qoiPath     = "Data/Temperature/Array"
computeMode = "THRESHOLD_ABS"
thresh      = 1100.0

IDT = AGATH_1D_IDT(solutDir, qoiPath, computeMode, threshold=thresh, abstol=1.0e-1, verbosity=1, plot=True)

prt("The End", "green", verbosity)
