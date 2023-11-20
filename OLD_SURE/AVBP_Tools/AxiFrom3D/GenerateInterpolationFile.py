import os
import numpy as np
import InterpolateFromAVBP as ci
import common_IO as cio
import math
from GenerateInterpolationFile_Axi    import *
from GenerateInterpolationFile_Cart   import *
from GenerateInterpolationFile_Slice  import *





def GenInterpFile(workDir, meshRelPath, interpRelPath, lim1, lim2, lim3, cart=False, Slice=False):
  if Slice:
      GenInterpFile_Slice(workDir, meshRelPath, interpRelPath, lim1, lim2)
  else:
    if cart:
      GenInterpFile_Cart(workDir, meshRelPath, interpRelPath, lim1, lim2, lim3)
    else:
      GenInterpFile_Axi(workDir, meshRelPath, interpRelPath, lim1, lim2, lim3)
