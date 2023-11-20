#!/usr/bin/env python2.7
import argparse
import Link
from AVBP_Tools.DecomposeMeans import DecomposeMeans


parser = argparse.ArgumentParser(description="Generates the decomposed means ready for bootstrapping")
parser.add_argument('-f', '--files', nargs='+', type=str , help = "List of solutions")
args = parser.parse_args()

files = args.files
if len(files) > 1:
    DecomposeMeans(files)

