import h5py as h5
import os
import numpy as np

def DecomposeMeans(lFile):
  nFile = len(lFile)
  lFile.sort()
  print(lFile[0] + " : keep as is")
  os.system("cp " + lFile[0] + " cut_" + lFile[0])
  for i in range(nFile-1):
    big = h5.File(lFile[i+1], 'r')
    lit = h5.File(lFile[i], 'r')
    tb  = big["Parameters/t_av"][()]
    tl  = lit["Parameters/t_av"][()]
    big.close()
    lit.close()
    if tb > tl:
      print(lFile[i+1] + " : decompose")
      _SingleDecomposeMeans(lFile[i+1], lFile[i])
    else:
      print(lFile[i+1] + " : keep as is")
      os.system("cp " + lFile[i+1] + " cut_" + lFile[i+1])


def _SingleDecomposeMeans(bigFile, littleFile):

  big       = h5.File(bigFile, 'r')
  lit       = h5.File(littleFile, 'r')
  key_names = lit["Average"].keys()
  tb        = big["Parameters/t_av"][()]
  tl        = lit["Parameters/t_av"][()]
  tn        = tb - tl

  new       = h5.File('cut_' + bigFile, 'w')
  new.create_group('Parameters')
  new.create_group('Average')

  # finding decomposed mean
  for key in key_names:
    fieldb = big['Average'][key][:]
    fieldl = lit['Average'][key][:]
    fieldn = np.zeros_like(fieldl)
    fieldn[:] = (tb * fieldb[:] - tl * fieldl[:]) / tn
    new.create_dataset('Average/'+key, data=fieldn)

  # filing parameters
  new.create_dataset('Parameters/t_av', data=tn)
  big.close()
  lit.close()
  new.close()
