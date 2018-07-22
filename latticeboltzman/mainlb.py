from latticeboltzman import *
import numpy as np 
import argparse as argparse;



parser = argparse.ArgumentParser(prog="python mainlb.py",
  description = "Initialises and configures a simple lattice boltzmann von karman sheet.");
parser.add_argument('--reynolds', '-R', help='Reynolds Number.', type=int, default=2500)

args = parser.parse_args();

reynolds = args.reynolds;

print "Chosen Re=%d" % reynolds;
simObj = lbsim(reynolds);
simObj.show ();