#!/usr/bin/env python2.7
import argparse
import h5py
import numpy as np




parser = argparse.ArgumentParser(description="Adds solved variance fields")
parser.add_argument('--file', '-f' , nargs='+', type=str , help = "List of files")
parser.add_argument('--var',  '-v' , nargs='+', type=str , help = "List of variables")
args = parser.parse_args()

files   = args.file
fields  = args.var

for filename in files:
  print(filename)
  f = h5py.File(filename,'r+')

  if "Zv_res" in fields:
    Z           = f['Average/YZ'][:]
    Z2          = f['Average/Y2Z'][:]
    newField    = np.zeros_like(Z)
    newField[:] = Z2[:] - (Z[:] * Z[:])
    try:
      f.create_dataset("Average/Zv_res", data=newField)
    except:
      del f["Average/Zv_res"]
      f.create_dataset("Average/Zv_res", data=newField)
    #fields.remove("Zv_res")
  #if "Ycv_res" in fields:
  #  Yc          = f['Average/YYc'][:]
  #  Yc2         = f['Average/Y2Yc'][:]
  #  newField    = np.zeros_like(Yc)
  #  newField[:] = Yc2[:] - (Yc[:] * Yc[:])
  #  try:
  #    f.create_dataset("Average/Ycv_res", data=newField)
  #  except:
  #    del f["Average/Ycv_res"]
  #    f.create_dataset("Average/Ycv_res", data=newField)
  #  #fields.remove("Ycv_res")
  #if "Chi_res" in fields:
  #  Chi         = f['Average/Chi_out'][:]
  #  Chi_sgs     = f['Average/Chi_sgs_out'][:]
  #  newField    = np.zeros_like(Chi)
  #  newField[:] = Chi[:] - Chi_sgs[:]
  #  try:
  #    f.create_dataset("Average/Chi_res", data=newField)
  #  except:
  #    del f["Average/Chi_res"]
  #    f.create_dataset("Average/Chi_res", data=newField)
  #  #fields.remove("Chi_res")
  #if "Zv_tot" in fields:
  #  try:
  #    Zv_res = f['Average/Zv_res'][:]
  #  except:
  #    Z           = f['Average/YZ'][:]
  #    Z2          = f['Average/Y2Z'][:]
  #    Zv_res      = np.zeros_like(Z)
  #    Zv_res[:]   = Z2[:] - (Z[:] * Z[:])
  #    f.create_dataset("Average/Zv_res", data=Zv_res)
  #  Zv_sgs      = f['Average/Yzv'][:]
  #  newField    = np.zeros_like(Chi)
  #  newField[:] = Zv_sgs[:] + Zv_res[:]
  #  try:
  #    f.create_dataset("Average/Zv_tot", data=newField)
  #  except:
  #    del f["Average/Zv_tot"]
  #    f.create_dataset("Average/Zv_tot", data=newField)
  #  #fields.remove("Zv_tot")

  #if "Z_Bilger" in fields:
  #  pass
  #  #YH          = f['Average/YH'][:]
  #  #YH          = f['Average/YH'][:]
  #  #Z2          = f['Average/Y2Z'][:]
  #  #Zv_res      = np.zeros_like(Z)
  #  #Zv_res[:]   = Z2[:] - (Z[:] * Z[:])
  #  #f.create_dataset("Average/Zv_res", data=Zv_res)



  #  #Zv_sgs      = f['Average/Yzv'][:]
  #  #newField    = np.zeros_like(Chi)
  #  #newField[:] = Zv_sgs[:] + Zv_res[:]
  #  #try:
  #  #  f.create_dataset("Average/Z_Bilger", data=newField)
  #  #except:
  #  #  del f["Average/Z_Bilger"]
  #  #  f.create_dataset("Average/Z_Bilger", data=newField)
  #  #fields.remove("Z_Bilger")

  f.close()


