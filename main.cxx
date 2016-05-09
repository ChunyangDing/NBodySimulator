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
#include <math.h>   // floor(), pow()
#include <iostream>
#include <fstream>
 #include <string.h> // memset()
#include <stdlib.h> // malloc, drand48()
#include <unistd.h> // getopt()
#include <cmath>    //fmod()

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
  double aMAX = .1;
  double da = 0.001;

  // Particle variables
  // double* x = new double[npart];
  // double* y = new double[npart];
  // double* z = new double[npart];

  double x[npart];
  double y[npart];
  double z[npart];

  // double *vx = new double[npart];
  // double *vy = new double[npart];
  // double *vz = new double[npart];

  double vx[npart];
  double vy[npart];
  double vz[npart];

  // double *gx = new double[ngrid * ngrid * ngrid];
  // double *gy = new double[ngrid * ngrid * ngrid];
  // double *gz = new double[ngrid * ngrid * ngrid];

  double gx[ngrid*ngrid*ngrid];
  double gy[ngrid*ngrid*ngrid];
  double gz[ngrid*ngrid*ngrid];

  double Dplus; // factor in initial conditions

  // real data has dimensions n*n*n
  double *rho = new double[ngrid * ngrid * ngrid];
  double *phi = new double[ngrid * ngrid * ngrid];

  // initialize rho and phi with 0s
  memset(rho, 0, sizeof(double)*ngrid*ngrid*ngrid);
  memset(phi, 0, sizeof(double)*ngrid*ngrid*ngrid);





  // setup initial conditions
  for (int i=0; i<npart; i++){
    // start with uniform distribution
    x[i]=drand48()*ngrid;
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

    vx[i]=0.5;  // initial velocities are all zero
    vy[i]=0.5;  
    vz[i]=0.5;  
  }
  
  //
  //////////////////////////////////////////////////////////////////////////
  while ( a < aMAX ){ 
    cicInterpolate(x, y, z, rho);

    solvePoisson(a, rho, phi);
    


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



    Field_on_Mesh(gx,gy,gz, phi); //solves for g on mesh


    updateParticles(a,da, x,y,z, vx,vy,vz, gx,gy,gz);


    a += da;
  }
  
    // delete[] x;
    // delete[] y;
    // delete[] z;

    // delete[] vx;
    // delete[] vy;
    // delete[] vz;

    // delete[] gx;
    // delete[] gy;
    // delete[] gz;

    delete[] rho;
    delete[] phi;

  outFile.close();


  
  return 0;
}
  
  


void cicInterpolate(double *x, double *y, double *z, double *rho){

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

  //for (int i=0; i<ngrid*ngrid*ngrid; i++){
  //  rho[i] = ( rho[i] - rhoCrit ) / rhoCrit;
  //}
}




//////////////////////////////////////////////////////////////////////////////////////////

void Field_on_Mesh(double *gx, double *gy, double *gz, double *phi)
{
  int ip1,im1,jp1,jm1,kp1,km1;
  int N=ngrid;

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
	
	gx[i*N*N + j*N +k]= -.5*(phi[ip1*N*N +  j*N  +  k ] - phi[im1*N*N + j*N  + k ]);
	gy[i*N*N + j*N +k]= -.5*(phi[ i*N*N  + jp1*N +  k ] - phi[ i*N*N + jm1*N + k ]);
	gz[i*N*N + j*N +k]= -.5*(phi[ i*N*N  +  j*N  + kp1] - phi[ i*N*N  + j*N + km1]); 

      }
    }
  } 
}

void solvePoisson(double a, double *rho, double *phi) {
    /* 1. declarations of relevant variables */ 
    int N=ngrid;
  fftw_complex *frho= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2 +1));
  memset(frho, 0, sizeof(fftw_complex)*N*N*(N/2+1));

  fftw_plan pf;
  fftw_plan pb;
  double kx, ky, kz;

  fftw_complex *fphi= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*(N/2 + 1));
  memset(fphi, 0, sizeof(fftw_complex)*N*N*(N/2+1));

  //double *green= new double[N*N*(N/2 +1)];
  double *green= (double*) malloc(sizeof(double)*N*N*(N/2 +1));
  memset(green, 0, sizeof(double)*N*N*(N/2+1));


    /* 2. fourier transform rho into fourier space
    We suggest looking to the fftw_plan_dft_r2c_3d() function for this transform. 
    Any fourier transform in fftw is then followed by the execution command. 
    For example:*/


    // setup fftw plan 
    pf = fftw_plan_dft_r2c_3d(ngrid, ngrid, ngrid, rho, frho, FFTW_ESTIMATE);

    // take fourier transform
    fftw_execute(pf); 


  /* 3. calculate green's function in fourier space */ 

    int xcounter=0;
    int ycounter=0;
    int zcounter=0;

    for (int l=-(N/2); l<(N/2); l++){
      ycounter = 0;
    for (int m=-(N/2); m<(N/2); m++){
      zcounter = 0;         
    for (int n=0; n<(N/2 +1); n++){
      //if ( (l==0) && (m==0) && (n==0)) green[xcounter*N*(N/2 +1) + ycounter*(N/2 +1) + zcounter-6] = 0;
      kx = 2.0*pi*l/N;
      ky = 2.0*pi*m/N;
      kz = 2.0*pi*n/N;
      // determine green
      if (xcounter*N*(N/2 +1) + ycounter*(N/2 +1) + zcounter < N*N*(N/2+1)) green[xcounter*N*(N/2 +1) + ycounter*(N/2 +1) + zcounter] = -((3.0*OmegaM)/(8.0*a))*pow( pow( sin(.5*kx),2) + pow(sin(.5*ky),2) + pow(sin(.5*kz),2),-1);
      else cout << xcounter << ' ' << ycounter << ' ' << zcounter << endl;
    zcounter++;}ycounter++;}xcounter++;}

    for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
    for (int k=0; k<(N/2 +1); k++){
      fphi[i*N*(N/2 +1) + j*(N/2 + 1) + k][0] = green[i*N*(N/2 +1) + j*(N/2 + 1) + k] * frho[i*N*(N/2 +1) + j*(N/2 + 1) + k][0];
      fphi[i*N*(N/2 +1) + j*(N/2 + 1) + k][1] = green[i*N*(N/2 +1) + j*(N/2 + 1) + k] * frho[i*N*(N/2 +1) + j*(N/2 + 1) + k][1];
    }}}


    /* 4. reverse transformation using fftw_plan_dft_c2r_3d() and fftw_execute() */
    pb = fftw_plan_dft_c2r_3d(ngrid, ngrid, ngrid, fphi, phi, FFTW_ESTIMATE);

    fftw_execute(pb);

    for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
    for (int k=0; k<N; k++){
      phi[i*N*N + j*N + k] *= 1.0/(N*N*N);
    }}} 

    //delete[] green;
    free(green);
    fftw_free(frho);
    fftw_free(fphi);
    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb);


    //(testing) brings back fftw to its initial state without any wisdom information retained
    //fftw_cleanup();
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

    cout << i << ' ' << j << ' ' << k << endl;

    int N = npart;

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
    vx[p] += ax*da;//f(a) * ax * da;
    vy[p] += ay*da;//f(a) * ay * da;
    vz[p] += az*da;//f(a) * az * da;
    /* update particle positions */ 
    x[p] += vx[p]*da + .5*ax*da*da;//pow((a+da/2.0), -2) * f(a + da/2.0) * vx[p] * da;
    y[p] += vy[p]*da + .5*ay*da*da;//pow((a+da/2.0), -2) * f(a + da/2.0) * vy[p] * da;
    z[p] += vz[p]*da + .5*az*da*da;//pow((a+da/2.0), -2) * f(a + da/2.0) * vz[p] * da;

    /* Handles wrap around particles, either to the left or to the right */
    // if (x[p] > ngrid) x[p] = x[p] - (double) N;
    // if (y[p] > ngrid) y[p] = y[p] - (double) N;
    // if (z[p] > ngrid) z[p] = z[p] - (double) N;
    // x[p] = fmod(x[p], ngrid);
    // y[p] = fmod(y[p], ngrid);
    // z[p] = fmod(z[p], ngrid);

    // if (x[p] < 0) x[p] = ngrid + x[p];
    // if (y[p] < 0) y[p] = ngrid + y[p];
    // if (z[p] < 0) z[p] = ngrid + z[p];
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
