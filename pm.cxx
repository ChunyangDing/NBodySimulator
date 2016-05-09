///////////////////////////////////////////////////////////////////////////////
/**  pm.cxx ::  PM N-body code
 *   
 *   Authors: Joseph, Chunyang, Mehmet
 *	 Date   : May 9. 2016
 *   *********************************
 *	
 *   This code models the gravitational interaction of many bodies in a
 *   collisionless system via the Particle-Mesh method.
 *
 *   FILES:
 *
 *   	pm.cxx is the main function.
 *   	pm.h   contains global variables and function prototypes
 *
 *   	data is written out to 'particle_positions.txt' and 'velocity.txt'
 *
 *   	This data is read by 'plotting.py' and 'plot_vx_vs_v.py' respectively
 *
 *   	These files generate plots of the data at each timestep which are 
 *   	saved in the directories 3D and VX.
 *
 *	 KNOWN-ISSUES:
 *
 *		D+ :: We don't have a clue what this is, but we have it set to 5
 *			  (in general the whole initial condition is a little shakey)
 *		
 *		Leap-frog Integration :: another thing we don't have a clue about
 *            (we believe we are currently using 1st-order Euler scheme)
 *
 *
 *
 */
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>    // pow(), sin()
#include <fftw3.h>
#include <cstdlib>  // srand48(), drand48()

#include "pm.h"

using namespace std;

int main()
{


// time-stepping variables
	  double a    = 0.010;
const double aMAX = 0.200;
const double da   = 0.001;

// particle arrays :: N = nparticle
double x[N];
double y[N];
double z[N];

double vx[N];
double vy[N];
double vz[N];

// grid arrays (fftw uses its own malloc())
double *rho;
double *phi;

rho = (double*) fftw_malloc(sizeof(double)*Mtot);
phi = (double*) fftw_malloc(sizeof(double)*Mtot);




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



/* setup initial conditions */


double Dplusdot=1; // basically useless, but don't change it
Dplus = 5;

/** 
 *  Particle locations are distributed randomly through out the 64x64x64 grid
 *  Then they are shifted in the Zel'dovich approximation. Only the x-component
 *  is shifted (1D-plane wave). We aren't totally confident about this.
 */


// set seed
srand48(1);

// loop over particles
for (int i=0; i<N; i++){
	x[i]=drand48()*M1D;
	y[i]=drand48()*M1D;
	z[i]=drand48()*M1D;



	x[i]= ((i-1)*M1D / N) + Dplus * sin(2*pi*((i-1)*M1D / N) / M1D);


	vx[i]=Dplusdot*sin(2*pi*x[i]/M1D);
	vy[i]=0;
	vz[i]=0;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/* set-up output files */

ofstream VXoutFile, out3D;
out3D.open("particle_positions.txt");
VXoutFile.open("velocity.txt");

/***********************/

// main

/* loop over time steps */
while ( a < aMAX )
{
	cicInterpolate(x,y,z, rho);

	solvePoisson(a, rho, phi);

	// print out data
	for (int i=0; i<N; i++){
		out3D   << x[i] << ' '
				<< y[i] << ' '
				<< z[i] << ' '
						<< '\n';
	}

	for (int i=0; i<N; i++){
		VXoutFile << x[i]  << ' '
				  << vx[i] << ' '
						   << '\n';
	}

	updateParticles(a,da, x,y,z, vx,vy,vz, phi);


	a += da;
}



/* close output files */
VXoutFile.close();
out3D.close();


/* free memory */
fftw_free(rho);
fftw_free(phi);


	return 0;
}










void cicInterpolate(double *x, double *y, double *z, double *rho){

  int I, J, K;      // Parent cell coordinates

  double dx,dy,dz;  // distances from parent cell center to particle

  double tx,ty,tz;  // interpolation factors


  /* loop over particles and apply the CIC interpolation to calculate density of each cell */
  for (int particle=0; particle<N; particle++)
  {
    I = floor(x[particle]);
    J = floor(y[particle]);
    K = floor(z[particle]);

    dx = x[particle] - (double) I;    tx = 1 - dx;
    dy = y[particle] - (double) J;    ty = 1 - dy;
    dz = z[particle] - (double) K;    tz = 1 - dz;


    rho[   (I %M1D)   *M1D*M1D +   (J %M1D)   *M1D +   (K %M1D)]   += tx*ty*tz;
    rho[   (I %M1D)   *M1D*M1D + ((J+1) %M1D) *M1D +   (K %M1D)]   += tx*dy*tz;
    rho[   (I %M1D)   *M1D*M1D +   (J %M1D)   *M1D + ((K+1) %M1D)] += tx*ty*dz;
    rho[   (I %M1D)   *M1D*M1D + ((J+1) %M1D) *M1D + ((K+1) %M1D)] += tx*dy*dz;
    rho[ ((I+1) %M1D) *M1D*M1D +   (J %M1D)   *M1D +   (K %M1D)]   += dx*ty*tz;
    rho[ ((I+1) %M1D) *M1D*M1D + ((J+1) %M1D) *M1D +   (K %M1D)]   += dx*dy*tz;
    rho[ ((I+1) %M1D) *M1D*M1D +   (J %M1D)   *M1D + ((K+1) %M1D)] += dx*ty*dz;
    rho[ ((I+1) %M1D) *M1D*M1D + ((J+1) %M1D) *M1D + ((K+1) %M1D)] += dx*dy*dz;
  }

}


void solvePoisson(double a, double *rho, double *phi) {
/**
 *   I am just going to say it. This function is a mess.
 *   I am 95% sure it is corrent.
 *
 *   The order that things are declared and executed is really important
 *   and the size of the arrays was a pain and a half to find out.
 *
 */





	fftw_complex *frho;
	const int frho_size = M1D*M1D*(M1D/2 + 1);
	frho = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*frho_size);

	/* fourier transform rho into fourier space */

	// create plan (let fftw3 learn machine architecture)
	fftw_plan pf;
	pf = fftw_plan_dft_r2c_3d(M1D, M1D, M1D, rho, frho, FFTW_ESTIMATE);

	fftw_execute(pf);
	fftw_destroy_plan(pf);



	/* calculate green's function in fourier space */
	double G[frho_size];
	double kx,ky,kz;

	for (int l=0; l<M1D; l++){
	for (int m=0; m<M1D; m++){
	for (int n=0; n<(M1D/2 + 1); n++){

		// frequencies
		kx = 2*pi*(l - (M1D/2)) / M1D;
		ky = 2*pi*(m - (M1D/2)) / M1D;
		kz = 2*pi*n / M1D;

		// done out here to save horizontal line space
		double sin2kx = pow(sin(kx/2.0), 2);
		double sin2ky = pow(sin(ky/2.0), 2);
		double sin2kz = pow(sin(kz/2.0), 2);

		// set (0,0,0) component to 0
		if ( (kx==0) && (ky==0) && (kz==0) ) {
			G[l*M1D*(M1D/2 + 1) + m*(M1D/2 + 1) + n] = 0;
		}
		else {
			G[l*M1D*(M1D/2 + 1) + m*(M1D/2 + 1) + n] = -((3.0*OmegaM)/(8.0*a))*pow( (sin2kx + sin2ky + sin2kz) ,-1);
		}
	}}}


	// get phi array ready

	fftw_complex *fphi;
	fphi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*frho_size);

	// make plan for backward transformation
	fftw_plan pb;
	pb = fftw_plan_dft_c2r_3d(M1D, M1D, M1D, fphi, phi, FFTW_ESTIMATE);

	// calculated necessary stuff
	for (int l=0; l<M1D; l++){
	for (int m=0; m<M1D; m++){
	for (int n=0; n<(M1D/2 + 1); n++){
			int pos=l*M1D*(M1D/2 + 1) + m*(M1D/2 + 1) + n;
		
		fphi[pos][0] = G[pos]*frho[pos][0];
		fphi[pos][1] = G[pos]*frho[pos][1];
			}
		}	
	}



	/* reverse transformation using fftw_plan_dft_c2r_3d() and fftw_execute() */

	fftw_execute(pb);
	fftw_destroy_plan(pb);

	/* normalize phi */
	for (int i=0; i<Mtot; i++){
		phi[i]= (1/(double)Mtot)*phi[i];
	}


	/* free memory */
	fftw_free(frho);
	fftw_free(fphi);
}



void updateParticles(double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *phi)
{
	/* declarations of relevant variables*/
	int ip1,im1,jp1,jm1,kp1,km1;

	double *gx;
	double *gy;
	double *gz;

	// potential ON THE MESH
	gx = (double*) malloc(sizeof(double)*Mtot);
	gy = (double*) malloc(sizeof(double)*Mtot);
	gz = (double*) malloc(sizeof(double)*Mtot);




	int I, J, K;      // Parent cell coordinates
	int IP1, JP1, KP1;

	double dx,dy,dz;  // distances from parent cell center to particle

	double tx,ty,tz;  // interpolation factors

	// accelerations of particle
	double ax;
	double ay;
	double az;


	// calculate gradient of potential ON THE MESH
	for (int i=0; i<M1D; i++){
	for (int j=0; j<M1D; j++){
	for (int k=0; k<M1D; k++){
		ip1 = i+1;  im1 = i-1;
		jp1 = j+1;  jm1 = j-1;
		kp1 = k+1;  km1 = k-1;

		// Boundary setting
		if (i == M1D-1) ip1 = 0;
		if (j == M1D-1) jp1 = 0;
		if (k == M1D-1) kp1 = 0;
		if (i == 0) im1 = M1D-1;
		if (j == 0) jm1 = M1D-1;
		if (k == 0) km1 = M1D-1;

		gx[i*M1D*M1D + j*M1D +k]= -.5*(phi[ip1*M1D*M1D +  j*M1D  +  k ] - phi[im1*M1D*M1D + j*M1D  + k ]);
		gy[i*M1D*M1D + j*M1D +k]= -.5*(phi[ i*M1D*M1D  + jp1*M1D +  k ] - phi[ i*M1D*M1D + jm1*M1D + k ]);
		gz[i*M1D*M1D + j*M1D +k]= -.5*(phi[ i*M1D*M1D  +  j*M1D  + kp1] - phi[ i*M1D*M1D  + j*M1D + km1]); 

	}}}


	/* calculate particle accelerations from grad phi */
	for (int p=0; p<N; p++){

		I = floor(x[p]);    IP1 = I+1;
		J = floor(y[p]);	JP1 = J+1;
		K = floor(z[p]);	KP1 = K+1;

		if (IP1 == M1D) IP1 = 0;
		if (JP1 == M1D) JP1 = 0;
		if (KP1 == M1D) KP1 = 0;

		dx = x[p] - (double) I;    tx = 1 - dx;
		dy = y[p] - (double) J;    ty = 1 - dy;
		dz = z[p] - (double) K;    tz = 1 - dz;

		ax = gx[I*M1D*M1D + J*M1D + K]*tx*ty*tz + gx[IP1*M1D*M1D + J*M1D + K]*dx*ty*tz +
			 gx[I*M1D*M1D + JP1*M1D + K]*tx*dy*tz + gx[IP1*M1D*M1D + JP1*M1D + K]*dx*dy*tz +
			 gx[I*M1D*M1D + J*M1D + KP1]*tx*ty*dz + gx[IP1*M1D*M1D + J*M1D + KP1]*dx*ty*dz +
			 gx[I*M1D*M1D + JP1*M1D + KP1]*tx*dy*dz + gx[IP1*M1D*M1D + JP1*M1D + KP1]*dx*dy*dz;

		ay = gy[I*M1D*M1D + J*M1D + K]*tx*ty*tz + gy[IP1*M1D*M1D + J*M1D + K]*dx*ty*tz +
			 gy[I*M1D*M1D + JP1*M1D + K]*tx*dy*tz + gy[IP1*M1D*M1D + JP1*M1D + K]*dx*dy*tz +
			 gy[I*M1D*M1D + J*M1D + KP1]*tx*ty*dz + gy[IP1*M1D*M1D + J*M1D + KP1]*dx*ty*dz +
			 gy[I*M1D*M1D + JP1*M1D + KP1]*tx*dy*dz + gy[IP1*M1D*M1D + JP1*M1D + KP1]*dx*dy*dz;

		az = gz[I*M1D*M1D + J*M1D + K]*tx*ty*tz + gz[IP1*M1D*M1D + J*M1D + K]*dx*ty*tz +
			 gz[I*M1D*M1D + JP1*M1D + K]*tx*dy*tz + gz[IP1*M1D*M1D + JP1*M1D + K]*dx*dy*tz +
			 gz[I*M1D*M1D + J*M1D + KP1]*tx*ty*dz + gz[IP1*M1D*M1D + J*M1D + KP1]*dx*ty*dz +
			 gz[I*M1D*M1D + JP1*M1D + KP1]*tx*dy*dz + gz[IP1*M1D*M1D + JP1*M1D + KP1]*dx*dy*dz;


		/* update particle velocities */
		vx[p] += f(a) * ax * da;
		vy[p] += f(a) * ay * da;
		vz[p] += f(a) * az * da;
		/* update particle positions */ 
		x[p] += pow((a+da/2.0), -2) * f(a + da/2.0) * vx[p] * da;
		y[p] += pow((a+da/2.0), -2) * f(a + da/2.0) * vy[p] * da;
		z[p] += pow((a+da/2.0), -2) * f(a + da/2.0) * vz[p] * da;

	}


	free(gx);
	free(gy);
	free(gz);


}

// helper function which prints out 3D arrays
void printVec3D(double* rho) {
  for (int i=0; i<M1D; i++) {
    for (int j=0; j<M1D; j++) {
      for (int k=0; k<M1D; k++) {
        cout << rho[i*M1D*M1D + j*M1D + k] << ' ';
        }
      cout << endl;
    }
  cout << endl;
  }
}

double f(double a){
  return pow((1 / a) * sqrt( OmegaM + OmegaK * a + OmegaL * pow(a, 2) ), -0.5);
}

