import matplotlib 
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.animation as animation 
from matplotlib.ticker import NullFormatter
import sys as sys;

class lbsim (object):
	"""simpleNavierStokes is an implementation of the ft08 fenics tutorial script.
	It encapsulates the script as an animation object, such that matplotlib displays
	the realtime calculations. In the future, this should just be done by inheriting
	some kind of animationObject superclass."""
	def __init__(self, reynoldsRequired):
		""" Initialise matplotlib functions, set some parameters""" 
 		self.animationInterval = 1; 
		self.figure = plt.figure( figsize=(20,20)); 
		self.time = 0;

		#sufficiently long length/time scales LB=NS

		#parameters
		self.snapshot   = 100;                    # Take a snapshot every nSnap time steps
		self.reynolds   = reynoldsRequired;      # Reynolds number (Re)

		scaler = 1.0;
		self.horizontalNodes  = int(520*scaler);      # Number of nodes in x direction
		self.verticalNodes    = int(180*scaler);      # Number of nodes in y direction
		self.channelWidth     = self.verticalNodes - 1.0; # Width of the channel

		self.horizontalCentre = self.horizontalNodes/4;     # Center of cylinder in x
		self.verticalCentre   = self.verticalNodes/2;     # Center of cylinder in y
		self.cylinderRadius   = self.verticalNodes/9;     # Cylinder radius

		self.velocityLatticeBoltzmann = 0.04;  # Velocity in lattice units
		# Kinematic viscosity associated with Re
		self.viscosityLatticeBoltzmann   = self.velocityLatticeBoltzmann*self.cylinderRadius/self.reynolds
		self.omega   = 1.0/(3.0*self.viscosityLatticeBoltzmann + 0.5)    # Relaxation parameter

		# 2d9Q velocity vectors
		self.velocityVectors = np.array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]);

		# Associated weigths: 4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36
		self.equilibriumWeights = (1.0/36.0)*np.ones(9);
		self.equilibriumWeights[np.asarray([np.linalg.norm(ci)<1.1 for ci in self.velocityVectors])] = 1.0/9.0;
		self.equilibriumWeights[0] = 4.0/9.0;

		# Bounce back boundary for no-slip condition on the cylinder
		self.noSlip = [self.velocityVectors.tolist().index((-self.velocityVectors[i]).tolist()) for i in range(9)] ;
 
		# Which of the velocities must be special cased
		# Unknown pops moving right
		self.specialCaseNegative = np.arange(9)[np.asarray([ci[0]<0  for ci in self.velocityVectors])];
		# Unknown pops in the middle
		self.specialCaseZero = np.arange(9)[np.asarray([ci[0]==0 for ci in self.velocityVectors])]; 
		# Unknown pops moving left
		self.specialCasePositive = np.arange(9)[np.asarray([ci[0]>0  for ci in self.velocityVectors])]; 


		# Mark the notes associated with the cylinder
		self.obstacle = np.fromfunction(lambda x,y: (x-self.horizontalCentre)**2+(y-self.verticalCentre)**2<self.cylinderRadius**2, (self.horizontalNodes, self.verticalNodes));

		# Impose a sinusoidal flow condition for the inlet and outlet
		# Actually, it is imposed everywhere in the first step
		self.initialVelocity = np.fromfunction(lambda d,x,y: (1-d)*self.velocityLatticeBoltzmann*(1.0 + 1.0e-4*np.sin(y/self.channelWidth*2*np.pi)),(2, self.horizontalNodes,self.verticalNodes));

		# Determine the initial equilibrium and 'previous step' values
		self.initialEquilibrium = self.equilibrium(1.0, self.initialVelocity);
		self.distributionsIn = self.initialEquilibrium.copy();
		self.distributionsOut = self.initialEquilibrium.copy();

	def sumPopulation(self, distributions):
		return np.sum(distributions, axis=0); 
	def equilibrium(self, density, velocities):

	    velocityComponents   = 3.0*np.dot(self.velocityVectors, velocities.transpose(1,0,2));

	    velocityLength = (3.0/2.0)*(velocities[0]**2 + velocities[1]**2)
	    equilibriumDistribution = np.zeros((9, self.horizontalNodes,self.verticalNodes))
	    for i in range(9):
			equilibriumDistribution[i,:,:] = density*self.equilibriumWeights[i]*(1.0 + velocityComponents[i] + 0.5*velocityComponents[i]**2 - velocityLength)
	    return equilibriumDistribution

	def square(self, p):
 
		return np.multiply(p,p);

 	def evolve(self, iteration):   

		# Right edge: the unknown populations are determined
		# from the known populations in the adjacent cell
		# resulting in an outflow boundary condition
		self.distributionsIn[self.specialCaseNegative,-1,:] = self.distributionsIn[self.specialCaseNegative,-2,:]; 

		# Calculate the density from the populations
		rho = self.sumPopulation(self.distributionsIn)

		# Calculate the velocity from the populations
		velocity = np.dot(self.velocityVectors.transpose(), self.distributionsIn.transpose((1,0,2)))/rho 

		velocity[:,self.obstacle] = 0.0 

		# Impose the original sinusoidal flow condition on the left-most cells 
		velocity[:,0,:] = self.initialVelocity[:,0,:]

		# Left edge: compute density from known populations
		rho[0,:] = 1./(1.-velocity[0,0,:]) * (self.sumPopulation(self.distributionsIn[self.specialCaseZero,0,:])+2.*self.sumPopulation(self.distributionsIn[self.specialCaseNegative,0,:]))

		# Update the equilibrium populations
		equilibriumDistribution = self.equilibrium(rho, velocity)

		# Left edge: adjust the populations to their imposed equilibrium values
		# this constitutes Zou/He boundary conditions (only the unknown
		# populations are modified, unlike in regular bounce-back, which modifies 
		# them all)
		self.distributionsIn[self.specialCasePositive,0,:] = equilibriumDistribution[self.specialCasePositive,0,:]

		# Collision step
		self.distributionsOut = self.distributionsIn - self.omega * (self.distributionsIn - equilibriumDistribution)

		# Impose the bounce-back no-slip condition
		for i in range(9): 
			self.distributionsOut[i, self.obstacle] = self.distributionsIn[self.noSlip[i],self.obstacle]

		# Stream the populations according to the LB velocities
		for i in range(9):
			self.distributionsIn[i,:,:] = np.roll(np.roll(self.distributionsOut[i,:,:], self.velocityVectors[i,0],axis=0), self.velocityVectors[i,1],axis=1)

		# Visualize the magnitude of the flow field
		self.time = self.time + 1;
		if ( self.time == 3 or self.time%self.snapshot == 0 ):
			plt.clf(); 
			plt.imshow(np.sqrt(velocity[0]**2 + velocity[1]**2).transpose(),cmap=cm.afmhot)
			plt.xlabel("x")
			plt.ylabel("y")
			plt.title("Karman Vortex Street t=%d reynolds=%d"% (self.time, self.reynolds)) 
 
	def show(self): 
		plt.plot([], [], 'r-')
		self.ani = animation.FuncAnimation(self.figure, self.evolve, interval=self.animationInterval) 
		plt.show ();
		