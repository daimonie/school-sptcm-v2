import numpy as np
import sys as sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm 
import argparse as argparse;
import time as time;

#pre-allocate random numbers (it is faster)
uniformContainer = np.random.rand(1000000); 
uniformCounter = 0;
def getRandom():
	global uniformCounter;
	global uniformContainer;

	uniformCounter = uniformCounter + 1;
	if uniformCounter > 1000000-1:
		uniformContainer = np.random.rand(1000000); 
		uniformCounter = 1;
	return uniformContainer[uniformCounter-1];


def gillespieFunc(number, length, startPopulation, mu):
 


	# Master equation
	#
	#  \partial_t P(N, t) = (N-1) gamma P(N-1, t) 
	#						+ mu P(N+2, t) (N+2) (N+1)/2
	#						- N gamma P(N,t)
	#						- mu P(N, t) (N)(N-1))/2
	#
	#	interested in statistics involving ratios of gamma/mu, so set gamma=1

	 
	populations = np.zeros((number, length));

	extinction = 0;
	for i in range(0, number): 
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
	 			#print "Population %d went extinct." % i; 
	 			extinction = extinction + 1;

	return extinction;

startTime = time.time();

sizeMu = 25;
sizeLength = 25;

muArray = np.linspace(0.01, .25, sizeMu);
lengthArray = np.linspace(0, 1e4, sizeLength);

 
final = sizeMu * sizeLength; 

progress = 0;

[muGrid, lengthGrid] = np.meshgrid(muArray, lengthArray);

extinction = 0 * muGrid;

for mm in range(0, sizeMu):
	for ll in range(0, sizeLength):

		length = lengthArray[ll];
		mu = muArray[mm]; 

		extinction[ll, mm] = gillespieFunc(1000, int(length), 1000, mu);

		print "Extinction %d, progress %d/%d, elapsed time %.3f" % ( extinction[ll, mm], progress, final, time.time() - startTime); 

		progress = progress + 1;

print "Elapsed time %.3f" % (time.time() - startTime);
# plot results
fig = plt.figure();
ax = fig.gca(projection='3d');
surf = ax.plot_surface(muGrid, lengthGrid, extinction, cmap=cm.afmhot);

ax.set_xlabel('mu');
ax.set_ylabel('length');
ax.set_zlabel('#Extinctions');

plt.show ();