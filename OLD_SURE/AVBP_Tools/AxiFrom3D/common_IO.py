import numpy as np
import os
import sys
import h5py as h5

def load_interp_base(filename):
  try:
    data = np.load(filename)
    print "Reading vertices and weights"
    return data
  except IOError:
    print "\nError!\n\n Dude, you've gotta provide a correct base file!"
    print " Check <fname_interp_base> variable or pathname!"
    print " Hint: The file is expected to be in your current working directory by default."
    sys.exit(1)


def save_interp_base(filename, vtx, wts, header, has_voln, **voln):
  try:
    if has_voln==True:
      np.savez_compressed(filename, vtx=vtx, wts=wts, voln=voln, header=header)
      print "Dumped vertices, weights and cell volumes to file for reuse!"
    else:
      np.savez_compressed(filename, vtx=vtx, wts=wts, header=header)
      print "Dumped vertices and weights to file for reuse!"
      print "You can add cell volumes for fancier images to your AVBP solution."
      print "Just add `save_solution.user_variable = voln` to your run.params"
      print "and run for 1 iteration."
  except IOError:
    print "\nError!\n\n"
    print "Can't write the file! Check your permissions or pathnames!"
    sys.exit(1)


def save_cutplane(filename, plane, header):
  try:
    np.savez_compressed(filename, plane=plane, header=header)
    print "Dumped cutplane for reuse!"
  except IOError:
    print "\nError!\n\n"
    print "Can't write the file! Check your permissions or pathnames!"
    sys.exit(1)



def save_sol(filename, field, header):
  try:
    np.savez_compressed(filename, field=field, header=header)
    print "Dumped solution for reuse!"
  except IOError:
    print "\nError!\n\n"
    print "Can't write the file! Check your permissions or pathnames!"
    sys.exit(1)

def load_sol(filename):
  try:
    data = np.load(filename)
    print "Reading vertices and weights"
    return data
  except IOError:
    print "\nError!\n\n Dude, you've gotta provide a correct base file!"
    print " Check <fname_interp_base> variable or pathname!"
    print " Hint: The file is expected to be in your current working directory by default."
    sys.exit(1)

def read_avbp_mesh_coordinates(filename):
  try:
    print "Reading mesh file"
    mesh = h5.File(filename,'r')
    
    print "Extracting coordinates from meshfile"
    x = mesh.get('Coordinates/x')
    x = np.array(x)
    y = mesh.get('Coordinates/y')
    y = np.array(y)
    z = mesh.get('Coordinates/z')
    z = np.array(z)
    mesh.close()
    return x, y, z
  except IOError:
    print "\nError!\n\n"
    print "Are you SERIAL, bro??? Better check your mesh path!"
    sys.exit(1)


def read_avbp_solution(filename, varname):
  try:
    f = h5.File(filename,'r')
    print "Reading solution file %s" % (filename.split('/')[-1])
    
    print "Extracting " + varname + " from solution file"
    old_var = varname
    if varname in f.keys():
      solution = f.get(varname)
    elif 'Average' in f.keys():
      if varname in f['Average'].keys():
        solution = f.get('Average/' + varname)
      else:
        solution = f.get(varname)
    elif 'Additionals' in f.keys():
      if varname in f['Additionals'].keys():
        solution = f.get('Additionals/' + varname)
      elif varname in f['FictiveSpecies'].keys():
        solution = f.get('FictiveSpecies/' + varname)
      elif varname in f['GaseousPhase'].keys():
        solution = f.get('GaseousPhase/' + varname)
      elif varname in f['RhoSpecies'].keys():
        solution = f.get('RhoSpecies/' + varname)
      else:
        solution = f.get(varname)
    else:
      solution = f.get(varname)
  
    solution = np.array(solution)
    if solution.shape == ():
      raise ValueError()	
    f.close()
    
    return solution
  except IOError:
    print "\nError!\n\n"
    print "You don't have this file, bro: " + filename
    print "Better check your path and filename!"
    sys.exit(1)
  except ValueError as e:
    print "\nError!\n\nCan't find variable " + varname + " in solution file " + filename
    print "Check your <VARKEY> entry and make sure the variable exists in the corresponding file."	
    sys.exit(1)


def write_avbp_solution(filename, varname, array):
  print("Opening solution file %s" % (filename.split('/')[-1]))
  f = h5.File(filename,'r+')
  try:
    data = f[varname][:]
    print("Overwriting field") # TODO: add yes/no
    f[varname][:] = array
  except KeyError:
    print("Creating dataset")
    f.create_dataset(varname, data=array)
  f.close()




def check_path_save(pathname):
  if os.path.exists(pathname):
    return True
  else:
    print "\nError!"
    print "You don't have this path, bro: %s" % pathname
    sys.exit(1)
