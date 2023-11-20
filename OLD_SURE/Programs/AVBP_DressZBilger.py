#!/usr/bin/env python2.7
import argparse
import Link
from AVBP_Tools.Dress.ZBilger import DressZBilgerCabra

parser = argparse.ArgumentParser(description="Dresses ZBilger in AVBP average sol")
parser.add_argument('-s', '--sol'   , type=str , help = "Solution file. Default is 'av_sol.h5'")
args = parser.parse_args()

if args.sol is not None:
  sol     = args.sol
else:
  sol     = 'av_sol.h5'

DressZBilgerCabra(sol)
