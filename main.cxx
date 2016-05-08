//////////////////////////////////////////////////////////////////////////////
/**
 *  main.cxx : Paricle Mesh N-body Code
 *
 *  Authors  : Chunyang, Joseph, Mehmet
 *  Date     : May 9, 2016
 *
 *  Uses PM technique to simulate gravitational interaction of N bodies in a
 *  collisionless system.
 *  ****************************
 *  Outline:
 *
 *    define location of particles according to Zel'dovich approximation
 *
 *    interpolate density at mesh points
 *
 *    determine gravitational potential at mesh points
 *
 *    interpolate gravitational potential at particle positions
 *
 *    calculate x,y,z forces on particle from potential
 *
 *    update particle momentum
 *
 *    update particle position
 *
 *
 *
 */
///////////////////////////////////////////////////////////////////////////////
//
//Current main error: The phi is off the walls insane. Huge numbers, and the first two rows are nonsense.
//Appears to be something wrong with the transform. Needs help with the Poisson portion. 
//


#include <fftw3.h>
#include <math.h> // floor(), pow()
#include <iostream>
#include <fstream>
#include <stdlib.h> // malloc, drand48()
#include <unistd.h> // getopt()
#include <cmath> //fmod()

#include "pm.h"

using namespace std;

int main(int argc, char* argv[])
{
  /////////////////////////////////////////////////////////////////////////////
  //File input and argument handling region
  int c; 
  while ((c = getopt(argc, argv,"Vh")) != -1){
    switch (c) {
      case 'V':
        cout << "Version of code: " << version << endl;
        return 0; //quit the program if 'V' is passed.
      case 'h':
        cout << "This program will simulate dark matter based on given input"
                "conditions and constants.  It will create a text file at the"
                "end with the coordinates and velocities of each particle at"
                "each time-step.  Please call the program in the following way:"
                 << endl;
        usage(argv[0]);
        cout << "All the arguments in the program are optional.\n-V\tprints"
                "the program version and quits\n-h\tprints instructions for"
                "the program and quits.\n[fileName]\t is an optional argument"
                "to specify the name of the output file."
                << endl;
        
        return 0; //quit the program if 'h' is passed.
      case '?':
        //call usage function if unkown option
        cout << "You have passed an unspecified option.  Run code as '"
             << argv[0] << " -h' for help." << endl;
				usage(argv[0]);
        return 2;
      }
  }
    
  //now check for the file name
  int optLeftOver = argc - optind;
  if (optLeftOver > 1) {
    cout << "You entered more than one standard argument or did not put [file]"
            " at the end.  Run code as '"
         << argv[0]
         << " -h' for help."
         << endl;
    return 2;
  }
    
  ofstream outFile;
  //if there is one opt left then that will be the name
  if (optLeftOver == 1) outFile.open(argv[optind]);
  else outFile.open("out.txt");

  //Make a header for the outputfile
  outFile << "ID\tx\ty\tz\tvx\tvy\tvz" << endl;


  /////////////////////////////////////////////////////////////////////////////
  //Creating local variables to be used throughout the code

  // time-stepping variables
  double a = 0.01;    // start when universe was 1% current size
  double aMAX = 0.1;  // End when universe is 100% current size
  double da = 0.01;   // Advance by 1% per time-step 

  //Giving me errors about "variable array length" even though npart is declared as constant.
  // Particle variables
  double* x = new double[npart];
  double* y = new double[npart];
  double* z = new double[npart];

  double *vx = new double[npart];
  double *vy = new double[npart];
  double *vz = new double[npart];
  
  double *gx = new double[ngrid * ngrid * ngrid];
  double *gy = new double[ngrid * ngrid * ngrid];
  double *gz = new double[ngrid * ngrid * ngrid];

  double Dplus; // factor in initial conditions

  double *myRho;
  double *myPhi;

  // real data has dimensions n*n*n
  myRho = (double*) fftw_malloc(sizeof(double)*ngrid*ngrid*ngrid);
  myPhi = (double*) fftw_malloc(sizeof(double)*ngrid*ngrid*ngrid);

  //DEPRECATED, FOR DELETION
  //printVec3D(ngrid, myPhi);
  //vector<double> rho;    // Should describe mass density for each cell
  //vector<double> phi;     // Should have unique density for each cell 
  
  fftw_complex *frho;
  fftw_complex *fphi;

  // complex arrays have dimension n*n*(n/2 + 1)
  frho= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ngrid*ngrid*(0.5*ngrid+1));
  fphi= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ngrid*ngrid*(0.5*ngrid+1));

  // setup initial conditions
  for (int i=0; i<npart; i++){
    // start with uniform distribution
    x[i]=drand48()*ngrid; //EDITED BY CHUNNY: should not go to L, but to ngrid. Bounds are 0 - 63, not 0 - 99
    y[i]=drand48()*ngrid;
    z[i]=drand48()*ngrid; 

    /* displace according to Zel'dovich
     *
     *  x = x0 + S(x0)
     *
     *  in our case (as per Darryl's suggestion)
     *  we have S(x) = sin( 2*pi*x/L )
     *
     *  D+ = 1 (although perhaps D+ = a)
     */
     Dplus = a;

    x[i] = x[i] + Dplus*sin( 2*pi*x[i] / L);
    y[i] = y[i] + Dplus*sin( 2*pi*y[i] / L);
    z[i] = z[i] + Dplus*sin( 2*pi*z[i] / L);

    vx[i]=0;  // initial velocities are all zero
    vy[i]=0;  
    vz[i]=0;  
  }
  
  //
  //////////////////////////////////////////////////////////////////////////
  //Main loop of the code; calculates the updates each step along the way
  // loop over time-steps
  while ( a < aMAX ){ 
    cout << a << endl;
    cicInterpolate(x, y, z, myRho);
    
    solvePoisson(a, myRho, frho, fphi, myPhi);
    
    outFile << '\n' << endl;
    for (int i=0; i<npart; i++) {
      outFile << '\n'
	      << x[i]  << "  \t"
	      << y[i]  << "  \t"
	      << z[i]  << "  \t"
	      << vx[i] << "  \t"
	      << vy[i] << "  \t"
	      << vz[i];
      }
    Field_on_Mesh(gx,gy,gx, myPhi); //solves for g on mesh
    updateParticles(a, da,x, y, z, vx, vy, vz, gx,gy,gz);
    a += da;
    //printVec3D(ngrid, myPhi);
  }
  
    
  outFile.close();
  fftw_free(frho);
  fftw_free(fphi);
  fftw_free(myRho);
  fftw_free(myPhi);
  
  return 0;
}
  
  

/*
 * Calculate contribution of each particle to the density grid points using the CIC interpolation.
 * Each particle contributes to the cell it is in and the neighboring grid cells based on the
 * algorithm in Section 2.8 of the write-up.
 */ 

void cicInterpolate(double *x, double *y, double *z, double *rho)
{
  /* declarations of relevant variables*/

  int I, J, K;      // Parent cell coordinates

  double dx,dy,dz;  // distances from parent cell center to particle

  double tx,ty,tz;  // interpolation factors


  /* loop over particles and apply the CIC interpolation to calculate density of each cell */
  for (int particle=0; particle<npart; particle++)
  {
    I = floor(x[particle]);
    J = floor(y[particle]);
    K = floor(z[particle]);

    dx = x[particle] - (double) I;    tx = 1 - dx;
    dy = y[particle] - (double) J;    ty = 1 - dy;
    dz = z[particle] - (double) K;    tz = 1 - dz;


    rho[   (I %ngrid)   *ngrid*ngrid +   (J %ngrid)   *ngrid +   (K %ngrid)]   += tx*ty*tz;
    rho[   (I %ngrid)   *ngrid*ngrid + ((J+1) %ngrid) *ngrid +   (K %ngrid)]   += tx*dy*tz;
    rho[   (I %ngrid)   *ngrid*ngrid +   (J %ngrid)   *ngrid + ((K+1) %ngrid)] += tx*ty*dz;
    rho[   (I %ngrid)   *ngrid*ngrid + ((J+1) %ngrid) *ngrid + ((K+1) %ngrid)] += tx*dy*dz;
    rho[ ((I+1) %ngrid) *ngrid*ngrid +   (J %ngrid)   *ngrid +   (K %ngrid)]   += dx*ty*tz;
    rho[ ((I+1) %ngrid) *ngrid*ngrid + ((J+1) %ngrid) *ngrid +   (K %ngrid)]   += dx*dy*tz;
    rho[ ((I+1) %ngrid) *ngrid*ngrid +   (J %ngrid)   *ngrid + ((K+1) %ngrid)] += dx*ty*dz;
    rho[ ((I+1) %ngrid) *ngrid*ngrid + ((J+1) %ngrid) *ngrid + ((K+1) %ngrid)] += dx*dy*dz;
  }


  // solve for delta = rho - rhoCrit / rhoCrit
  //********* This part is unclear *************//

  for (int i=0; i<ngrid*ngrid*ngrid; i++){
    rho[i] = ( rho[i] - rhoCrit ) / rhoCrit;
  }


}


/*  Solve poisson's equations and calculate acceleration field  */ 
void solvePoisson(double a, double *myRho, fftw_complex *frho,fftw_complex *fphi, double *myPhi)
{
  /* The FFT's are NOT done "in-place" (in the future we probably should switch)
   * Separate arrays are used for the input and the FFT'd data.
   *
   * The inverse (c2r) transform overwrites input fphi (doesn't matter)
   */
  /* declarations of relevant variables */
  double kx,ky,kz; // frequencies
  double *G = new double[ngrid * ngrid * (int) floor((0.5 * ngrid + 1))]; //Green's Function
  fftw_plan p_rho;
  fftw_plan p_phi;

  /* Lets do something stupid... the fftw3 library is meant for C, not C++ which
   * means that it doesn't want to be given vecotrs...
   * So I have to copy over the data from the vectors
   * into new arrays (and later copy them back)
   */

  //This is commented out because we no longer need to use vectors, and always use 1D double arrays
  // double *myRho;
  // double *myPhi;

  // // real data has dimensions n*n*n
  // myRho = (double*) fftw_malloc(sizeof(double)*ngrid*ngrid*ngrid);
  // myPhi = (double*) fftw_malloc(sizeof(double)*ngrid*ngrid*ngrid);
  // since we malloc'ed we have to free the memory ourselves later

  // for (int i=0; i<ngrid*ngrid*ngrid; i++){
  //   myRho[i] = rho[i];
  //   myPhi[i] = phi[i];
  // }


  // says to go from rho to frho
  p_rho = fftw_plan_dft_r2c_3d(ngrid,ngrid,ngrid, myRho, frho, FFTW_ESTIMATE);

  // take fourier transform (actually does it)
  fftw_execute(p_rho);
  fftw_destroy_plan(p_rho); // free some memory

  /* calculate green's function in fourier space */
  for (int l=0; l<ngrid; l++){
    for (int m=0; m<ngrid; m++){
      for (int n=0; n<(0.5*ngrid+1); n++){
	if ( (l==0) && (m==0) && (n==0) ) {G[0] = 0;}
	else {
	  kx = pi*(l-0.5*ngrid) / L;
	  ky = pi*(m-0.5*ngrid) / L;
	  kz = pi*n / L;
	  G[l+m+n] = -((3.0*OmegaM)/(8.0*a)) * pow( (pow(sin(kx),2) + pow(sin(ky),2) + pow(sin(kz),2)), -1);
	}
      }
    }
  }

  /* calculate fphi in fourier space */
  for (int l=0; l<ngrid; l++){
    for (int m=0; m<ngrid; m++){
      for (int n=0; n<0.5*ngrid+1; n++){
	fphi[l+m+n][0] = G[l+m+n]*frho[l+m+n][0];
	fphi[l+m+n][1] = G[l+m+n]*frho[l+m+n][1];
      }
    }
  }


  /* reverse transformation */

  p_phi = fftw_plan_dft_c2r_3d(ngrid,ngrid,ngrid, fphi, myPhi, FFTW_ESTIMATE);

  fftw_execute(p_phi);
  fftw_destroy_plan(p_phi);

  /* nomralize phi and copy back into vector*/
  for (int i=0; i<ngrid; i++){
    for (int j=0; j<ngrid; j++){
      for (int k=0; k<ngrid; k++){
	myPhi[i+j+k]= (1.0/pow(L,3))*myPhi[i+j+k];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void Field_on_Mesh(double *gx, double *gy, double *gz, double *phi)
{
  int ip1,im1,jp1,jm1,kp1,km1;

  for (int i=0; i<ngrid; i++){
    for (int j=0; j<ngrid; j++){
      for (int k=0; k<ngrid; k++){
	ip1 = i+1;  im1 = i-1;
	jp1 = j+1;  jm1 = j-1;
	kp1 = k+1;  km1 = k-1;

	// Boundary setting
	if (i == ngrid-1) ip1 = 0;
	if (j == ngrid-1) jp1 = 0;
	if (k == ngrid-1) kp1 = 0;
	if (i == 0) im1 = ngrid-1;
	if (j == 0) jm1 = ngrid-1;
	if (k == 0) km1 = ngrid-1;
	
	gx[i*ngrid*ngrid+j*ngrid+k]= -.5*(phi[ip1*ngrid*ngrid + j*ngrid + k] - phi[im1*ngrid*ngrid + j*ngrid + k]);
	gy[i*ngrid*ngrid+j*ngrid+k]= -.5*(phi[i*ngrid*ngrid + jp1*ngrid + k] - phi[i*ngrid*ngrid + jm1*ngrid + k]);
	gz[i*ngrid*ngrid+j*ngrid+k]= -.5*(phi[i*ngrid*ngrid + j*ngrid + kp1] - phi[i*ngrid*ngrid + j*ngrid + km1]); 
	
      }
    }
  } 
}



void updateParticles(double a, double da, double *x, double *y, double*z, double *vx, double *vy, double*vz, double *gx, double *gy, double *gz){
  for (int p=0; p<npart; p++)
  {
    int i = static_cast<int>(floor(x[p]));
    int j = static_cast<int>(floor(y[p]));
    int k = static_cast<int>(floor(z[p]));

    double dx = x[p] - i;
    double dy = y[p] - j;
    double dz = z[p] - k;

    double tx = 1 - dx;
    double ty = 1 - dy; 
    double tz = 1 - dz;

    double ax;
    double ay;
    double az;
    // Boundary conditions
    int ip1, jp1, kp1;
    ip1 = i+1;
    jp1 = j+1;
    kp1 = k+1;

    if (i == ngrid-1) ip1 = 0;
    if (j == ngrid-1) jp1 = 0;
    if (k == ngrid-1) kp1 = 0;

    int N = ngrid;
    //cout << i << ", " << j << ", " << k << endl;
    ax= gx[i*N*N +  j*N  + k ]*tx*ty*tz + gx[ip1*N*N +  j*N + k ]*dx*ty*tz + 
      gx[i*N*N + jp1*N + k ]*tx*dy*tz + gx[ip1*N*N + jp1*N+ k ]*dx*dy*tz + 
      gx[i*N*N +  j*N + kp1]*tx*ty*dz + gx[ip1*N*N +  j*N +kp1]*dx*ty*dz + 
      gx[i*N*N + jp1*N +kp1]*tx*dy*dz + gx[ip1*N*N+ jp1*N + kp1]*dx*dy*dz;
    
    ay= gy[i*N*N +  j*N  + k ]*tx*ty*tz + gy[ip1*N*N +  j*N + k ]*dx*ty*tz + 
      gy[i*N*N + jp1*N + k ]*tx*dy*tz + gy[ip1*N*N + jp1*N+ k ]*dx*dy*tz + 
      gy[i*N*N +  j*N + kp1]*tx*ty*dz + gy[ip1*N*N +  j*N +kp1]*dx*ty*dz + 
      gy[i*N*N + jp1*N +kp1]*tx*dy*dz + gy[ip1*N*N+ jp1*N + kp1]*dx*dy*dz;

    az= gz[i*N*N +  j*N  + k ]*tx*ty*tz + gz[ip1*N*N +  j*N + k ]*dx*ty*tz + 
      gz[i*N*N + jp1*N + k ]*tx*dy*tz + gz[ip1*N*N + jp1*N+ k ]*dx*dy*tz + 
      gz[i*N*N +  j*N + kp1]*tx*ty*dz + gz[ip1*N*N +  j*N +kp1]*dx*ty*dz + 
      gz[i*N*N + jp1*N +kp1]*tx*dy*dz + gz[ip1*N*N+ jp1*N + kp1]*dx*dy*dz;

    /* update particle velocities */
    vx[p] += f(a) * ax * da;
    vy[p] += f(a) * ay * da;
    vz[p] += f(a) * az * da;
    /* 4. update particle positions */ 
    x[p] += pow((a+da/2.0), -2) * f(a + da/2.0) * vx[p] * da;
    y[p] += pow((a+da/2.0), -2) * f(a + da/2.0) * vy[p] * da;
    z[p] += pow((a+da/2.0), -2) * f(a + da/2.0) * vz[p] * da;

    /* Handles wrap around particles, either to the left or to the right */
    x[p] = fmod(x[p], ngrid);
    y[p] = fmod(y[p], ngrid);
    z[p] = fmod(z[p], ngrid);

    if (x[p] < 0) x[p] = ngrid + x[p];
    if (y[p] < 0) y[p] = ngrid + y[p];
    if (z[p] < 0) z[p] = ngrid + z[p];
  }
}

















///////////////////////////////////////////////////////////////////////////////////////////

// //Some helper functions for calculating the g for the x, y, and z directions 
// double getGx(int i, int j, int k, double *phi){
//   double gx;
//   if (i == 0){
//     gx = -0.5 * (phi[ngrid*ngrid* (ngrid - 1) + ngrid * j + k] - phi[ngrid*ngrid * (i + 1) + ngrid * j + k]);// i = ngrid - 1
//   }
//   else {
//     if (i == ngrid){
//       gx = - 0.5 * (phi[ngrid*ngrid * (i - 1) + ngrid * j + k] - phi[ngrid * j + k]); //i = 0  
//     }
//     else{
//       gx = -0.5 * (phi[ngrid*ngrid * (i - 1) + ngrid * j + k] - phi[ngrid*ngrid * (i + 1) + ngrid * j + k]); // The normal state
//     }
//   }
//   return gx;
// }


// double getGy(int i, int j, int k, double *phi){
//   double gy;
//   if (j == 0){
//     gy = -0.5 *(phi[ngrid*ngrid* i + ngrid * (ngrid - 1) + k] - phi[ngrid*ngrid * i + ngrid * (j+1) + k]);// j = ngrid - 1
//   }
//   else {
//     if (j == ngrid){
//       gy = -0.5 * (phi[ngrid*ngrid * i + ngrid * (j - 1) + k] - phi[ngrid*ngrid * i + k]); // j = 0
//     }
//     else{
//       gy = -0.5 * (phi[ngrid*ngrid * i + ngrid * (j - 1) + k] - phi[ngrid*ngrid * i + ngrid * (j + 1) + k]); // normal
//     }
//   }
//   return gy;
// }



// double getGz(int i, int j, int k, double *phi){
//   double gz;
//   if (k == 0){
//     gz = -0.5 * (phi[ngrid*ngrid * i + ngrid * j + (ngrid - 1)] - phi[ngrid*ngrid * + ngrid * j + (k + 1)]); // k = ngrid - 1
//   }
//   else {
//     if (k == ngrid){
//       gz = -0.5 * (phi[ngrid*ngrid * i + ngrid * j + (k - 1)] - phi[ngrid*ngrid * i + ngrid * j]); // k = 0
//     }
//     else {
//       gz = -0.5 * (phi[ngrid*ngrid * i + ngrid * j + (k - 1)] - phi[ngrid*ngrid * i + ngrid * j + (k + 1)] ); //normal
//     }
//   }
//   return gz;
// }




// /* Update position, velocities for each particle */ 
// void updateParticles(double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *phi)
// {
//   for (int abc = 0; abc < npart; abc++){
//     /* 1. declarations of relevant variables */ 
//     int i = static_cast<int>(floor(*x));
//     int j = static_cast<int>(floor(*y));
//     int k = static_cast<int>(floor(*z));

//     double dx = *x - i;
//     double dy = *y - i;
//     double dz = *z - i;

//     double tx = 1 - dx;
//     double ty = 1 - dy; 
//     double tz = 1 - dz;

//     double gx = 0;
//     double gy = 0;
//     double gz = 0;

//     /* 2. calculate particle accelerations from phi */
//     gx = getGx(i, j, k, phi) * tx * ty * tz + getGx(i+1, j, k, phi) * dx * ty * tz + getGx(i, j+1, k, phi) * tx * dy * tz + getGx(i + 1, j+1, k, phi) * dx * dy * tz + getGx(i, j, k+1, phi)* tx * ty * dz + getGx(i+1, j, k+1, phi) * dx * ty * dz + getGx(i, j+1, k+1, phi) * tx * dy * dz + getGx(i + 1, j+1, k+1, phi) * dx * dy * dz;
//     gy = getGy(i, j, k, phi) * tx * ty * tz + getGy(i+1, j, k, phi) * dx * ty * tz + getGy(i, j+1, k, phi) * tx * dy * tz + getGy(i + 1, j+1, k, phi) * dx * dy * tz + getGy(i, j, k+1, phi)* tx * ty * dz + getGy(i+1, j, k+1, phi) * dx * ty * dz + getGy(i, j+1, k+1, phi) * tx * dy * dz + getGy(i + 1, j+1, k+1, phi) * dx * dy * dz;
//     gz = getGz(i, j, k, phi) * tx * ty * tz + getGz(i+1, j, k, phi) * dx * ty * tz + getGz(i, j+1, k, phi) * tx * dy * tz + getGz(i + 1, j+1, k, phi) * dx * dy * tz + getGz(i, j, k+1, phi)* tx * ty * dz + getGz(i+1, j, k+1, phi) * dx * ty * dz + getGz(i, j+1, k+1, phi) * tx * dy * dz + getGz(i + 1, j+1, k+1, phi) * dx * dy * dz;

//     /* 3. update particle velocities */
//     vx[abc] = vx[abc] + f(a) * gx * da;
//     vy[abc] = vy[abc] + f(a) * gy * da;
//     vz[abc] = vz[abc] + f(a) * gz * da;
//     /* 4. update particle positions */ 
//     x[abc] = x[abc] + pow((a+da/2.0), -2) * f(a + da/2.0) * vx[abc] * da;
//     y[abc] = y[abc] + pow((a+da/2.0), -2) * f(a + da/2.0) * vy[abc] * da;
//     z[abc] = z[abc] + pow((a+da/2.0), -2) * f(a + da/2.0) * vz[abc] * da;


//     while (x[abc] > 100) x[abc] -= L;
//     while (y[abc] > 100) y[abc] -= L;
//     while (z[abc] > 100) z[abc] -= L;
//   }
// }





double f(double a){
  return pow((1 / a) * sqrt( OmegaM + OmegaK * a + OmegaL * pow(a, 2) ), -0.5);
}

void usage(const char* prog) {
    cerr << "Usage: " << prog << " [-V] [-h] [fileName]" << endl;
}

void printVec3D(int ngrid, double* rho) {
  for (int i=0; i<ngrid; i++) {
    for (int j=0; j<ngrid; j++) {
      for (int k=0; k<ngrid; k++) {
        cout << rho[i*ngrid*ngrid + j*ngrid + k] << ' ';
        }
      cout << endl;
    }
  cout << endl;
  }
}
