# SURE
# Software for Uncertainty Range Evaluation
A UQ code based on OpenTURNS and scikit-learn


WARNING: Make sure you add the environment variable $SURE_HOME to your path before trying to launch a program.


To be able to organize all the python files in folders, every folder has to be turned into a package by the addition of the __init__.py in each of them.
Furthermore, so that all the modules can call each other, the programs must begin with :
import link
Which effectively does :
import sys
sys.path.insert(0, <path to $SURE_HOME>)
This adds the main "package" to the path in which python will search for modules afterwards


Development was done with python 2.7

Dependencies are (with versions used for development in parenthesis):
h5py         (2.9.0)
numpy        (1.15.4)
scipy        (1.2.0)
matplotlib   (2.2.3)
openturns    (1.13)
scikit-learn (0.20.3)


optional:
termcolor (1.1.0)
progress  (1.5) (With pip only)


Thanks to Typhaine Lavabre for the awesome name 
