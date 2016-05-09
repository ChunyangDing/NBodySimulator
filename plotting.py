###############################################################################
#   plotting.py :: Plot the particles in 3D
#
#	Reads in from file 'particle_positions.txt'
#	Does some fancy stuff to put the values in arrays, then makes plots for
#	each timestep, saved in a folder called 3D
#
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ntimestep=189 #hard-coded
npart = 32*32 #hard-coded


x = np.ndarray((npart))
y = np.ndarray((npart))
z = np.ndarray((npart))

# open file
with open("particle_positions.txt", 'r') as f:
	for t in range(ntimestep):
		for i in range(npart):
			line = f.readline().split()
			x[i] = float(line[0])
			y[i] = float(line[1])
			z[i] = float(line[2])

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		ax.set_xlim(0, 70)
		ax.set_ylim(0, 70)
		ax.set_zlim(0, 70)

		ax.scatter(x, y, z)


		plt.savefig("./3D/{0}.png".format(t))
		plt.close()






