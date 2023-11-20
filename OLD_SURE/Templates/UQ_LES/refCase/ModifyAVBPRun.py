import h5py
import time
import numpy as np
import sys
import os
import subprocess

import h5py


###################
##    Parsing     #
###################
#parser = argparse.ArgumentParser(description="Finish UFPV tabulation for AVBP")
#parser.add_argument('-t', '--table', required=True, help = "Path to table")
#parser.add_argument('--plot'   , action='store_true', help = "plot Frond")
#parser.add_argument('--noWrite', action='store_true', help = "Prevent effective encryption on table file")
#parser.add_argument('--dummySR', required=True, help = "SR at which the flame cannot ignite. A dummy flamelet is created at this SR with same values except source terms which are set to 0")
#args = parser.parse_args()
#
#table     = args.table
#plot      = args.plot
#noWrite   = args.noWrite
#dummySR   = float(args.dummySR)

workDir = "."

tableName = "UFPV_AVBP.h5"
table = workDir + "/TABLE/" + tableName
"""
Co-flow temperature can easily be read in the new UFPV table : it is the first value of the temperature array ! \o/
"""
t    = h5py.File(table, 'r')
cofT  = int(t["Data/TEMPERATURE"][0])
t.close()

###############################
#    Changing target table    #
###############################
# Modifying input file
mixtdat = workDir + "/RUN/mixture_database.dat"
txtfile = open (mixtdat, "r")  
lines = txtfile.readlines()
txtfile.close()
lines[7] = "  lookup_table = ../TABLE/" + tableName + "\n" # Table path is given on the 8th line of the mixture database file
with open(mixtdat, "w") as txtfile:
    for line in lines:
        txtfile.write(line)

###############################
#     Recomputing thermo      #
###############################
cmdList = ["mpirun", "-n", "1", "change_table_ttc"]
p = subprocess.Popen(cmdList, cwd=workDir + '/RUN/SOLUT')
p.wait()
print("toto !")

###############################
#    Modifying boundaries     #
###############################
# Modifying input file
gsbchoice = workDir + "/SOLUTBOUND/gensolutbound.choices"
txtfile = open (gsbchoice, "r")  
lines = txtfile.readlines()
txtfile.close()
lines[59] = str(cofT) + "d0\n" # Co flow temperature is specified on the 60th line of the gensolutbound file
with open(gsbchoice, "w") as txtfile:
    for line in lines:
        txtfile.write(line)

# Executing gensolutbound
cmdList = ["mpirun","-n", "1", "gensolutbound", "--step2"]
p = subprocess.Popen(cmdList, cwd=workDir + '/SOLUTBOUND')
p.wait()
# Resetting turbulence profile
p = subprocess.Popen(["python", "ModifyBoundary.py"], cwd=workDir + '/CUSTOMINLET')
p.wait()
print("The End !")
