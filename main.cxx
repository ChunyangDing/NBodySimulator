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
 */
///////////////////////////////////////////////////////////////////////////////

#include <fftw3.h>
#include <math.h> // floor(), pow()
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h> // malloc, drand48()
#include <unistd.h> // getopt()

#include "pm.h"

using namespace std;

int main(int argc, char* argv[])
{
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
  //if there is one opt left then that will be the name, else the file will be named by default
  if (optLeftOver == 1) outFile.open(argv[optind]);
  else outFile.open("out.txt");

  //Make a header for the outputfile
  outFile << "ID\tx\ty\tz\tvx\tvy\tvz" << endl;


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////



  // time-stepping variables
  double a = 0.01;    // start when universe was 1% current size
  double aMAX = 1.0;  // End when universe is 100% current size
  double da = 0.01;   // Advance by 1% per time-step 

  // Particle variables
  double x[npart];
  double y[npart];
  double z[npart];

  double vx[npart];
  double vy[npart];
  double vz[npart];
  



  double *myRho;
  double *myPhi;

  // real data has dimensions n*n*n
  myRho = (double*) fftw_malloc(sizeof(double)*ngrid*ngrid*ngrid);

  //vector<double> rho;    // Should describe mass density for each cell
  vector<double> phi;     // Should have unique density for each cell 

  fftw_complex *frho;
  fftw_complex *fphi;

  // complex arrays have dimension n*n*(n/2 + 1)
  frho= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ngrid*ngrid*(0.5*ngrid+1));
  fphi= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ngrid*ngrid*(0.5*ngrid+1));
    






  // setup initial conditions
  for (int i=0; i<npart; i++){
    // start with uniform distribution
    x[i]=drand48()*L;
    y[i]=drand48()*L;
    z[i]=drand48()*L; 

    /* displace according to Zel'dovich
     *
     *  x = x0 + S(x0)
     *
     *  in our case (as per Darryl's suggestion)
     *  we have S(x) = sin( 2*pi*x/L )
     *
     *  D+ = 1 (although perhaps D+ = a)
     */
    x[i] = x[i] + 10*sin( 2*pi*x[i] / L);
    y[i] = y[i] + 10*sin( 2*pi*y[i] / L);
    z[i] = z[i] + 10*sin( 2*pi*z[i] / L);

    vx[i]=0;  // initial velocities are all zero
    vy[i]=0;  
    vz[i]=0;  
  }
 
      outFile << "\nAt a = " << a;
      for (int i=0; i<npart; i++) {
        outFile << '\n'
                << x[i]  << " "
                << y[i]  << " "
                << z[i]  << " "
                << vx[i] << " "
                << vy[i] << " "
                << vz[i];
      }


  
  // // loop over time-steps
  // while ( a < aMAX )
  //   { 
      cicInterpolate(x, y, z, myRho);
      
      // solvePoisson(a, rho, frho, fphi, phi);
      
      // outFile << "\nAt a = " << a;
      // for (int i=0; i<npart; i++) {
      //   outFile << '\n'
      //           << x[i]  << " "
      //           << y[i]  << " "
      //           << z[i]  << " "
      //           << vx[i] << " "
      //           << vy[i] << " "
      //           << vz[i];
      // }
      
      // updateParticles(a, da,&x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], phi);
      
      // a += da;
    // }


  outFile.close();
  fftw_free(frho);
  fftw_free(fphi);
  fftw_free(myRho);

    return 0;
}



/*
 * Calculate contribution of each particle to the density grid points using the CIC interpolation.
 * Each particle contributes to the cell it is in and the neighboring grid cells based on the
 * algorithm in Section 2.8 of the write-up.
 */ 
// void cicInterpolate(double *x, double *y, double *z, vector<double>& rho) {
//     //declare parent cell locations
//     int pcx, pcy, pcz;
    
//     //declare distances from particle to cell (note that t_x/y/z is not explicitly defined but only used as '1 - d_x/y/x')
//     double dx, dy, dz;
    
//     //declare the factors that will be used for mass assignment (either d or t)
//     double xfactor, yfactor, zfactor;
    
//     //loop over particles
//     for (int counter=0; counter<npart; counter++) {
//         //get parent cell locations
//         pcx = floor(x[counter]);
//         pcy = floor(y[counter]);
//         pcz = floor(z[counter]);
        
//         dx = x[counter] - pcx;
//         dy = y[counter] - pcy;
//         dz = z[counter] - pcz;
        
//         //loop over relevant cells
//         for (int i=pcx; i<pcx+2; i++) {
//             //get the d and t variable for x
//             //if on the first iteration of x set the xfactor to tx, on the second iteration to dx
//             if (i - pcx == 0) {xfactor = 1 - dx;}
//             else {xfactor = dx;}
//             for (int j=pcy; j<pcy+2; j++) {
//                 //get the d and t variable for y
//                 //if on the first iteration of y set the yfactor to ty, on the second iteration to dy
//                 if (j - pcy == 0) {yfactor = 1 - dy;}
//                 else {yfactor = dy;}
//                 for (int k=pcz; k<pcz+2; k++) {
//                     //get the d and t variable for z
//                     //if on the first iteration of z set the zfactor to tz, on the second iteration to dz
//                     if (k - pcz == 0) {zfactor = 1 - dz;}
//                     else {zfactor = dz;}
                    
//                     //now increment the mass of the appropriate rho index (assuming mass = 1)
//                     rho[(i%ngrid)*ngrid*ngrid + (j%ngrid)*ngrid + (k%ngrid)] += xfactor*yfactor*zfactor;
//                 }
//             }
//         }
//     }
// }
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

}


/*  Solve poisson's equations and calculate acceleration field  */ 
void solvePoisson(double a, vector<double> &rho, fftw_complex *frho,fftw_complex *fphi, vector<double> &phi)
{
  /* The FFT's are NOT done "in-place" (in the future we probably should switch)
   * Separate arrays are used for the input and the FFT'd data.
   *
   * The inverse (c2r) transform overwrites input fphi (doesn't matter)
   */
  /* declarations of relevant variables */
  double kx,ky,kz; // frequencies
  double G[ngrid*ngrid*(int)floor((0.5*ngrid +1))]; // Green's function
  fftw_plan p_rho;
  fftw_plan p_phi;

  /* Lets do something stupid... the fftw3 library is meant for C, not C++ which
   * means that it doesn't want to be given vecotrs...
   * So I have to copy over the data from the vectors
   * into new arrays (and later copy them back)
   */
  double *myRho;
  double *myPhi;

  // real data has dimensions n*n*n
  myRho = (double*) fftw_malloc(sizeof(double)*ngrid*ngrid*ngrid);
  myPhi = (double*) fftw_malloc(sizeof(double)*ngrid*ngrid*ngrid);
  // since we malloc'ed we have to free the memory ourselves later

  for (int i=0; i<ngrid*ngrid*ngrid; i++){
    myRho[i] = rho[i];
    myPhi[i] = phi[i];
  }


  // says to go from rho to frho
  p_rho = fftw_plan_dft_r2c_3d(ngrid,ngrid,ngrid, myRho, frho, FFTW_ESTIMATE);

  // take fourier transform (actually does it)
  fftw_execute(p_rho);
  fftw_destroy_plan(p_rho); // free some memory

  /* calculate green's function in fourier space */
  for (int l=0; l<ngrid; l++){
  for (int m=0; m<ngrid; m++){
  for (int n=0; n<0.5*ngrid+1; n++){
    if ( (l==0) && (m==0) && (n==0) ) {G[0] = 0;}
    else {
      kx = pi*(l-0.5*ngrid) / L;
      ky = pi*(m-0.5*ngrid) / L;
      kz = pi*n / L;
  G[l+m+n] = -((3.0*OmegaM)/(8.0*a)) * pow( (pow(sin(kx),2) + pow(sin(ky),2) + pow(sin(kz),2)), -1);
    }
  }}}

  /* calculate fphi in fourier space */
  for (int l=0; l<ngrid+1; l++){
  for (int m=0; m<ngrid+1; m++){
  for (int n=0; n<0.5*ngrid+1; n++){
    fphi[l+m+n][0] = G[l+m+n]*frho[l+m+n][0];
    fphi[l+m+n][1] = G[l+m+n]*frho[l+m+n][1];
  }}}


  /* reverse transformation */

  p_phi = fftw_plan_dft_c2r_3d(ngrid,ngrid,ngrid, fphi, myPhi, FFTW_ESTIMATE);

  fftw_execute(p_phi);
  fftw_destroy_plan(p_phi);

  /* nomralize phi and copy back into vector*/
  for (int i=0; i<ngrid; i++){
  for (int j=0; j<ngrid; j++){
  for (int k=0; k<ngrid; k++){
    phi[i+j+k]= (1.0/pow(L,3))*myPhi[i+j+k];
  }}}


  fftw_free(myRho); fftw_free(myPhi);

}



//Some helper functions for calculating the g for the x, y, and z directions 
double getGx(int i, int j, int k, vector<double> &phi){
  double gx;
  if (i == 0){
    gx = -0.5 * (phi[pow(ngrid, 2)* (ngrid - 1) + ngrid * j + k] - phi[pow(ngrid, 2) * (i + 1) + ngrid * j + k]);// i = ngrid - 1
  }
  else {
    if (i == ngrid){
      gx = - 0.5 * (phi[pow(ngrid, 2) * (i - 1) + ngrid * j + k] - phi[ngrid * j + k]); //i = 0  
    }
    else{
      gx = -0.5 * (phi[pow(ngrid, 2) * (i - 1) + ngrid * j + k] - phi[pow(ngrid, 2) * (i + 1) + ngrid * j + k]); // The normal state
    }
  }
  return gx;
}


double getGy(int i, int j, int k, vector<double> &phi){
  double gy;
  if (j == 0){
    gy = -0.5 *(phi[pow(ngrid, 2)* i + ngrid * (ngrid - 1) + k] - phi[pow(ngrid, 2) * i + ngrid * (j+1) + k]);// j = ngrid - 1
  }
  else {
    if (j == ngrid){
      gy = -0.5 * (phi[pow(ngrid, 2) * i + ngrid * (j - 1) + k] - phi[pow(ngrid, 2) * i + k]); // j = 0
    }
    else{
      gy = -0.5 * (phi[pow(ngrid, 2) * i + ngrid * (j - 1) + k] - phi[pow(ngrid, 2) * i + ngrid * (j + 1) + k]); // normal
    }
  }
  return gy;
}



double getGz(int i, int j, int k, vector<double> &phi){
  double gz;
  if (k == 0){
    gz = -0.5 * (phi[pow(ngrid, 2) * i + ngrid * j + (ngrid - 1)] - phi[pow(ngrid, 2) * + ngrid * j + (k + 1)]); // k = ngrid - 1
  }
  else {
    if (k == ngrid){
      gz = -0.5 * (phi[pow(ngrid, 2) * i + ngrid * j + (k - 1)] - phi[pow(ngrid, 2) * i + ngrid * j]); // k = 0
    }
    else {
      gz = -0.5 * (phi[pow(ngrid, 2) * i + ngrid * j + (k - 1)] - phi[pow(ngrid, 2) * i + ngrid * j + (k + 1)] ); //normal
    }
  }
  return gz;
}




/* Update position, velocities for each particle */ 
void updateParticles(double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, vector<double> &phi)
{
  for (int abc = 0; abc < npart; abc++){
    /* 1. declarations of relevant variables */ 
    int i = static_cast<int>(floor(*x));
    int j = static_cast<int>(floor(*y));
    int k = static_cast<int>(floor(*z));

    double dx = *x - i;
    double dy = *y - i;
    double dz = *z - i;

    double tx = 1 - dx;
    double ty = 1 - dy; 
    double tz = 1 - dz;

    double gx = 0;
    double gy = 0;
    double gz = 0;

    /* 2. calculate particle accelerations from phi */
    gx = getGx(i, j, k, phi) * tx * ty * tz + getGx(i+1, j, k, phi) * dx * ty * tz + getGx(i, j+1, k, phi) * tx * dy * tz + getGx(i + 1, j+1, k, phi) * dx * dy * tz + getGx(i, j, k+1, phi)* tx * ty * dz + getGx(i+1, j, k+1, phi) * dx * ty * dz + getGx(i, j+1, k+1, phi) * tx * dy * dz + getGx(i + 1, j+1, k+1, phi) * dx * dy * dz;
    gy = getGy(i, j, k, phi) * tx * ty * tz + getGy(i+1, j, k, phi) * dx * ty * tz + getGy(i, j+1, k, phi) * tx * dy * tz + getGy(i + 1, j+1, k, phi) * dx * dy * tz + getGy(i, j, k+1, phi)* tx * ty * dz + getGy(i+1, j, k+1, phi) * dx * ty * dz + getGy(i, j+1, k+1, phi) * tx * dy * dz + getGy(i + 1, j+1, k+1, phi) * dx * dy * dz;
    gz = getGz(i, j, k, phi) * tx * ty * tz + getGz(i+1, j, k, phi) * dx * ty * tz + getGz(i, j+1, k, phi) * tx * dy * tz + getGz(i + 1, j+1, k, phi) * dx * dy * tz + getGz(i, j, k+1, phi)* tx * ty * dz + getGz(i+1, j, k+1, phi) * dx * ty * dz + getGz(i, j+1, k+1, phi) * tx * dy * dz + getGz(i + 1, j+1, k+1, phi) * dx * dy * dz;

    /* 3. update particle velocities */
    *vx = *(vx) + f(a) * gx * da;
    *vy = *(vy) + f(a) * gy * da;
    *vz = *(vz) + f(a) * gz * da;
    /* 4. update particle positions */ 
    *x = *x + pow((a+da/2.0), -2) * f(a + da/2.0) * *vx * da;
    *y = *y + pow((a+da/2.0), -2) * f(a + da/2.0) * *vy * da;
    *z = *z + pow((a+da/2.0), -2) * f(a + da/2.0) * *vz * da;

    //Move to the next particle
    x++;
    y++;
    z++;
    vx++;
    vy++;
    vz++;
  }
}





double f(double a){
  return pow((1 / a) * sqrt( OmegaM + OmegaK * a + OmegaL * pow(a, 2) ), -0.5);
}

void usage(const char* prog) {
    cerr << "Usage: " << prog << " [-V] [-h] [fileName]" << endl;
}
