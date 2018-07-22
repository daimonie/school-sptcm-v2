import numpy as np
import sys as sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse as argparse;
import time as time;

startTime = time.time();


parser = argparse.ArgumentParser(prog="python gillespie.py",
  description = "Simple gillespie stochastic simulation .");

parser.add_argument('--number', '-N', help='Number of calculations', type=int, default=10);
parser.add_argument('--length', '-L', help='Length of calculations', type=int, default=1000); 
parser.add_argument('--mu', '-M', help='Decay coefficient mu (gamma = 1)', type=float, default=0.01);
parser.add_argument('--start', '-S', help='Starting Population', type=int, default=1000);

args = parser.parse_args();
 
number = args.number;
length = args.length; 
startPopulation = args.start;
mu = args.mu;

#pre-allocate random numbers (it is faster)
uniformContainer = np.random.rand((number*length*2)); 
uniformCounter = 0;


# Master equation
#
#  \partial_t P(N, t) = (N-1) gamma P(N-1, t) 
#						+ mu P(N+2, t) (N+2) (N+1)/2
#						- N gamma P(N,t)
#						- mu P(N, t) (N)(N-1))/2
#
#	interested in statistics involving ratios of gamma/mu, so set gamma=1

def getRandom():
	global uniformCounter;

	uniformCounter = uniformCounter + 1;
	return uniformContainer[uniformCounter-1];
 
populations = np.zeros((number, length));

extinction = 0;
for i in range(0, number):
	if number%(1+number/100)==0:
		print "Simulating population %d/%d." % (i, number);
	for j in range(0, length):

		currentPopulation = startPopulation;
		if j > 0:
 			currentPopulation = populations[i, j-1];
		if currentPopulation <= 0:
			break;


 		rateProliferateUp = (currentPopulation - 1);
 		rateProliferateDown = currentPopulation;

 		rateDecayUp = mu * (currentPopulation+2)*(currentPopulation+1)  / 2.0;
 		rateDecayDown = mu * currentPopulation * (currentPopulation-1) / 2.0;

 		rateTotal = rateProliferateDown + rateProliferateUp + rateDecayDown + rateDecayUp;

 		dt = -1/rateTotal * np.log(1-getRandom ());

 		# "die" cast for which rate wins out
 		rateDie = rateTotal * getRandom ();

 		if rateDie < rateProliferateUp:
 			populations[i, j] = currentPopulation + 1;
 		elif rateDie < rateProliferateUp + rateProliferateDown:
 			populations[i, j] = currentPopulation + 1;
 		elif rateDie < rateTotal - rateDecayDown:
 			populations[i, j] = currentPopulation - 2;
 		else:
 			populations[i, j] = currentPopulation - 2;

 		if populations[i,j] <= 0:
 			print "Population %d went extinct." % i;
 			extinction = extinction + 1;

print " %d Populations went extinct, rate is %.3e" % (extinction, extinction/(1.0*number));

print "Elapsed time %.3f seconds." % (time.time() - startTime);
if number < 100:
	for i in range(number):
		plt.semilogy(populations[i, :]);
	plt.show (); 
