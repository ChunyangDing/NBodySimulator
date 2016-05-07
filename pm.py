from numpy import *
from math import floor
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


# function prototypes would go here


def main():

	# cosmological constants
	OmegaM  = 0.27;
	OmegaK  = 0.73;
	OmegaL  = 0.00008;
	H0      = 67.80;   		  #(km/s)/Mpc
	G       = 6.67*(10**-11); # N*m**2 / kg**2

	da = 0.1;  # time-stepping resolution
	a  = 0.01; # start when universe was 0.01 current size
	aMAX = 1.0 # end at current time (when universe is equal to current size)

	# mesh constants
	M = 64**1; # per direction
	L = 1.0;   # actual length of mesh in "suitable" units

	# particle constants
	N  = 32**1;
	mp = 1.0;

	# handy units
	r0 = L / M;    # size of grid cell
	t0 = H0**-1;   # unit of time
	v0 = r0 / t0;  # unit of velocity
	p0 = ( (3*H0**2)/(8*pi*G) ) * OmegaM; # unit of density
	phi0 = v0**2;  # unit of potential

	# Particle arrays
	x = empty(N);
	y = empty(N);
	z = empty(N);

	x_tilde = empty(N);   # don't technically need separate arrays for this
	y_tilde = empty(N);   # tilde arrays are for code units
	z_tilde = empty(N);   # doing this helps clarify the process

	vx = empty(N);
	vy = empty(N);
	vz = empty(N);

	px_tilde = empty(N);
	py_tilde = empty(N);
	pz_tilde = empty(N);

	# Grid arrays
	density       = zeros((M,M,M));
	density_tilde = zeros((M,M,M));
	delta_tilde   = zeros((M,M,M));
	delta_cap     = zeros((M,M,M));   # Fourier space
	Green         = zeros((M,M,M));   # Fourier space
	potential_cap = zeros((M,M,M));   # Fourier space
	potential_tilde = zeros((M,M,M));

	gx_tilde = zeros((M,M,M));   # at cell centers
	gy_tilde = zeros((M,M,M));
	gz_tilde = zeros((M,M,M));

	gx_p = zeros(N);   # acceleration at particle
	gy_p = zeros(N);
	gz_p = zeros(N);

	# assign particles
	for i in range(N):
		x[i] = random.uniform(L)/r0;	# randomly spaced
		y[i] = random.uniform(L)/r0;
		z[i] = random.uniform(L)/r0;

		vx[i] = 0; 				# no initial velocity
		vy[i] = 0;
		vz[i] = 0;

	###########################################################
	###########################################################
	# Plot initial

	fig = plt.figure();       # create figure object
	ax = fig.add_subplot(111, projection='3d'); # create 3D axes

	ax.set_xlim(-10,70);
	ax.set_ylim(-10,70);
	ax.set_zlim(-10,70);

	ax.scatter(x,y,z, c='g');  # plot particles

	ax.set_xlabel('X');
	ax.set_ylabel('Y');
	ax.set_zlabel('Z');

	plt.savefig("0.0.png");
	plt.close();
	############################################################
	############################################################

	t=0;
	TMAX=2;
	da = .01;
	# start while loop (for time)
	while t < TMAX:
		# compute density
		for particle in range(N):
			i = floor(x[particle]);
			j = floor(y[particle]);
			k = floor(z[particle]);

			dx = x[particle] - i;
			dy = y[particle] - j;
			dz = z[particle] - k;

			tx = 1 - dx;
			ty = 1 - dy;
			tz = 1 - dz;
			# Boundary conditions
			ip1 = i+1;
			jp1 = j+1;
			kp1 = k+1;

			if (i == M-1):
				ip1 = 0;
			if (j == M-1):
				jp1 = 0;
			if (k == M-1):
				kp1 = 0;

			density[ i,   j,   k ] += mp*tx*ty*tz;  	# these aren't densities...
			density[ i,  jp1,  k ] += mp*tx*dy*tz;
			density[ i,   j,  kp1] += mp*tx*ty*dz;
			density[ i,  jp1, kp1] += mp*tx*dy*dz;
			density[ip1,  j,   k ] += mp*dx*ty*tz;
			density[ip1, jp1,  k ] += mp*dx*dy*tz;
			density[ip1,  j,  kp1] += mp*dx*ty*dz;
			density[ip1, jp1, kp1] += mp*dx*dy*dz;

		# make density units right
		for i in range(M):
			for j in range(M):
				for k in range(M):
					density_tilde[i,j,k] = density[i,j,k] / p0;

		delta_tilde = density_tilde - 1;

		# Calculate Fourier Transformed density on the mesh
		delta_cap = fft.fftn(delta_tilde);

		# Calculate Green's function
		k = 2*pi / M;
		for l in range(M):
			for m in range(M):
				for n in range(M):
					if (l==0) and (m==0) and (n==0):
						Green[l,m,n] = 0;
					else:
						Green[l,m,n] = -( (3*OmegaM)/(8) )*( sin(.5*k*l)**2 + sin(.5*k*m)**2 + sin(.5*k*n)**2 )**-1;

		# Calculate (FFTed) potential on the mesh using Poisson's eqn
		for l in range(M):
			for m in range(M):
				for n in range(M):
					potential_cap[l,m,n] = Green[l,m,n]*delta_cap[l,m,n];


		# Find potential using an Inverse Fourier Transform
		potential_tilde = fft.ifftn(potential_cap);

		# for all particles do

			# compute field (acceleration) in nearby meshpoints
		for i in range(M):
			for j in range(M):
				for k in range(M):
					ip1 = i+1;  im1 = i-1;
					jp1 = j+1;  jm1 = j-1;
					kp1 = k+1;  km1 = k-1;

					# Boundary setting
					if (i == M-1):
						ip1 = 0;
					if (j == M-1):
						jp1 = 0;
					if (k == M-1):
						kp1 = 0;
					if (i == 0):
						im1 = M-1;
					if (j == 0):
						jm1 = M-1;
					if (k == 0):
						km1 = M-1;

					gx_tilde[i,j,k] = -( potential_tilde[ip1,j,k] - potential_tilde[im1,j,k] )/2;
					gy_tilde[i,j,k] = -( potential_tilde[i,jp1,k] - potential_tilde[i,jm1,k] )/2;
					gz_tilde[i,j,k] = -( potential_tilde[i,j,kp1] - potential_tilde[i,j,km1] )/2;

			# interpolate the field at the particle position
		for particle in range(N):
			i = floor(x[particle]);
			j = floor(y[particle]);
			k = floor(z[particle]);

			dx = x[particle] - i;
			dy = y[particle] - j;
			dz = z[particle] - k;

			tx = 1 - dx;
			ty = 1 - dy;
			tz = 1 - dz;

			# Boundary conditions
			ip1 = i+1;
			jp1 = j+1;
			kp1 = k+1;

			if (i == M-1):
				ip1 = 0;
			if (j == M-1):
				jp1 = 0;
			if (k == M-1):
				kp1 = 0;

			gx_p[particle] =  gx_tilde[i, j,  k ]*tx*ty*tz + gx_tilde[ip1, j,  k ]*dx*ty*tz + \
							  gx_tilde[i,jp1, k ]*tx*dy*tz + gx_tilde[ip1,jp1, k ]*dx*dy*tz + \
							  gx_tilde[i, j, kp1]*tx*ty*dz + gx_tilde[ip1, j, kp1]*dx*ty*dz + \
							  gx_tilde[i,jp1,kp1]*tx*dy*dz + gx_tilde[ip1,jp1,kp1]*dx*dy*dz;

			gy_p[particle] =  gy_tilde[i, j,  k ]*tx*ty*tz + gy_tilde[ip1, j,  k ]*dx*ty*tz + \
							  gy_tilde[i,jp1, k ]*tx*dy*tz + gy_tilde[ip1,jp1, k ]*dx*dy*tz + \
							  gy_tilde[i, j, kp1]*tx*ty*dz + gy_tilde[ip1, j, kp1]*dx*ty*dz + \
							  gy_tilde[i,jp1,kp1]*tx*dy*dz + gy_tilde[ip1,jp1,kp1]*dx*dy*dz;

			gz_p[particle] =  gz_tilde[i, j,  k ]*tx*ty*tz + gz_tilde[ip1, j,  k ]*dx*ty*tz + \
							  gz_tilde[i,jp1, k ]*tx*dy*tz + gz_tilde[ip1,jp1, k ]*dx*dy*tz + \
							  gz_tilde[i, j, kp1]*tx*ty*dz + gz_tilde[ip1, j, kp1]*dx*ty*dz + \
							  gz_tilde[i,jp1,kp1]*tx*dy*dz + gz_tilde[ip1,jp1,kp1]*dx*dy*dz;

			f = 1
			# calculate new speed and position using a suitable scheme
			vx[particle] = vx[particle] + f*gx_p[particle]*da;
			vy[particle] = vy[particle] + f*gy_p[particle]*da;
			vz[particle] = vz[particle] + f*gz_p[particle]*da;

			x[particle] = (x[particle] + f*vx[particle]*da) % 63;
			y[particle] = (y[particle] + f*vy[particle]*da) % 63;
			z[particle] = (z[particle] + f*vz[particle]*da) % 63;

		t = t + da;

		# plotting
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		ax.set_xlim(-10,70)
		ax.set_ylim(-10,70)
		ax.set_zlim(-10,70)

		ax.scatter(x,y,z, c='g')  # particles

		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		ax.set_zlabel('Z')

		plt.savefig("{0}.png".format(t))
		plt.close()

	return 0;


def f(a):
	return ( (a**-1) * (OmegaM + OmegaK*a + OmegaL*a**3) )**(-0.5);



if __name__ == __main__:
	main()