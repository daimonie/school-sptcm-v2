import numpy as np
import sys as sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse as argparse;



parser = argparse.ArgumentParser(prog="python langevin.py",
  description = "Simple Langevin trajectories .");

parser.add_argument('--start', '-S', help='Start Position Y.', type=float, default=0.0);
parser.add_argument('--number', '-N', help='Number of trajectories', type=int, default=100);
parser.add_argument('--length', '-L', help='Length of trajectories', type=int, default=1000);
parser.add_argument('--gamma', '-G', help='friction coefficient', type=float, default=0.5);

args = parser.parse_args();

start = args.start;
numberTrajectories = args.number;
lengthTrajectories = args.length;
gamma = args.gamma;




startPosition = np.array([0.0, start]);

mass = 1.0; 

forceParams = [1.0, 0.0, 0.0];
fluctuationParams = [0.0, 1.00, 10.0];
deltaTime = 0.01;

omega = 2 * np.pi / .1;
distance = lambda xx, yy: np.sqrt(xx*xx + yy*yy) + 1e-4;
forceL = lambda xx, rr: 2 * xx * np.cos(omega*rr) - rr * rr * np.sin(omega*rr) * xx / rr;

force = lambda xx, yy: forceParams[0]*np.array([forceL(xx, distance(xx,yy)), forceL(yy,distance(xx,yy))]) + np.array(forceParams[2:3]);
fluctuation = lambda mm, ss: fluctuationParams[2] * np.array([ np.random.normal(mm, ss, 1) for i in range(2)]);

trajectories = np.zeros((numberTrajectories, lengthTrajectories, 2));
velocities = np.zeros((numberTrajectories, lengthTrajectories, 2));

#as naive as possible
for i in range(0, numberTrajectories):
	for j in range(0, lengthTrajectories): 

		currentPosition = startPosition;
		currentVelocity = np.array([0.0, 0.0]);
		if j > 0:
			currentPosition = trajectories[i, j-1]; 
			currentVelocity = velocities[i, j-1]; 

		currentForce = force( currentPosition[0], currentPosition[1]);
		currentFluctuation = fluctuation(fluctuationParams[0], fluctuationParams[1]); 
 		#reshape
 		currentFluctuation = currentFluctuation.T[0];

		velocities[i,j] = (currentFluctuation + currentForce)/mass*deltaTime + (1.0 - gamma  * deltaTime) * currentVelocity;
		trajectories[i,j] = velocities[i, j] * deltaTime + currentPosition; 
 
fig = plt.figure();
ax = fig.add_subplot(111, projection='3d');

z = np.array( range(0, lengthTrajectories));

for i in range(0, numberTrajectories):
	ax.plot(trajectories[i, :, 0], trajectories[i, :, 1], z)



ax.set_xlabel('X position')
ax.set_ylabel('Y position')
ax.set_zlabel('Nth time step')

ax.set_title('%d trajectories of length %d '  % (numberTrajectories, lengthTrajectories));


xlimiter = np.max(np.max( np.abs(trajectories[i, :,0])))*10;

ax.set_xlim(-xlimiter, xlimiter);

ax.set_ylim(-xlimiter, xlimiter);
if start > xlimiter:
	ax.set_ylim(-1, start+1.0);

plt.show();