import matplotlib 
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation 
from matplotlib.ticker import NullFormatter

class wormlikechain (object):
	"""simpleNavierStokes is an implementation of the ft08 fenics tutorial script.
	It encapsulates the script as an animation object, such that matplotlib displays
	the realtime calculations. In the future, this should just be done by inheriting
	some kind of animationObject superclass."""
	def __init__(self):
		""" Initialise matplotlib functions, set some parameters"""
 
 		self.animationInterval = 1; 
		self.figure = plt.figure( figsize=(20,20));
		self.polymer = np.array([[0,0]]);
		self.length = 0;

	def square(self, p):
 
		return np.multiply(p,p);

 	def evolve(self, iteration):  

 		accepted = False;
 		newBead = np.array([0,0]);

 		while accepted == False:
	 		angle = np.random.rand() * 2 * np.pi; 

	 		distance = np.random.normal(0,1,1);  

	 		diff = distance*np.array([np.cos(angle), np.sin(angle)]); 

	 		lastBead = self.polymer[self.length];

	 		newBead = lastBead + diff; 

	 		accepted = np.abs(newBead[1]) < 3;

 		self.polymer = np.concatenate((self.polymer, [newBead]), axis=0);
 
 		self.length = self.length + 1;
 		if self.length%100==0:
			plt.plot(self.polymer[:,0], self.polymer[:,1]);

	 		Rsq = self.square(self.polymer[:, 0]) + self.square(self.polymer[:, 1]); 
	 		
	 		R = np.sqrt(np.max(Rsq));

	 		print "%d %.3f %.3f" % (self.length, np.sqrt(self.length), R);
 
	def show(self): 
		plt.plot([], [], 'r-')
		self.ani = animation.FuncAnimation(self.figure, self.evolve, interval=self.animationInterval) 
		plt.show ();
		