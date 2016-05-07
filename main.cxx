#include <fftw3.h>
#include <math.h> // for floor(), pow()
#include <iostream>
#include <vector>
#include <stdlib.h>

using namespace std;

//Don't forget to create function templates up here!
void cicInterpolate(int ngrid, int npart, double *x, double *y, double *z, double ***rho);
void solvePoisson(double a, double ***rho, double ***phi, int ngrid);
void updateParticles(int ngrid, int npart, double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, double ***phi);
int f(double a);

int main() 
{
  // relevant variables
  const int ngrid = pow(64,3);
  const int npart = pow(32,3);

  // cosmological variables
  const double OmegaM = 0.27;
  const double OmegaL = 0.73;
  const double OmegaK = 0.0000824;
  const double G = 6.6740831 * pow(10, -11);
  const double rhoCrit = 1.06 * pow(10, -29);
  const double H = 67.80;
  const double pi = 3.141592654;



  double x[npart]; //Should describe every particle
  double y[npart];
  double z[npart];
  
  double vx[npart];
  double vy[npart];
  double vz[npart];
  


  // variables for converting to actual units (we can ignore these for our purposes)
  //const double r0 = 1;
  //const double t0 = 1/H;
  //const double v0 = r0 / t0;
  //const double rho0 = ((3 * pow(H, 2)) / (8 * pi * G)) * OmegaM;
  //const double phi0 = pow(v0, 2);

  //Should be vectors now!
  //double rho[ngrid][ngrid][ngrid];    // Should describe mass density for each cell
  //double phi[ngrid][ngrid][ngrid];     // Should have unique density for each cell
  

  // int *** phi;
  // phi = new int**[ngrid];
  // for (int i = 0; i < ngrid; i++){
  //   phi[i] = new int*[ngrid];
  //   for (int j = 0; j<ngrid; j++){
  //     phi[i][j] = new int[ngrid];
  //   }
  // }

  // int *** rho;
  // rho = new int**[ngrid];
  // for (int i = 0; i < ngrid; i++){
  //   rho[i] = new int*[ngrid];
  //   for (int j = 0; j<ngrid; j++){
  //     rho[i][j] = new int[ngrid];
  //   }
  // }

  
  double a  = 0.1;     // Should be the scale parameter, usually = 1
  double aMAX = 1.0
  double da = 0.01;    // Some small change correlated with time step.
  
    
  /* 2. read in or setup initial conditions */
  
  for (int i=0; i<ngrid; i++)
    {
      x[i]=1;
      y[i]=1;
      z[i]=1;
      
      vx[i]=0;     
      vy[i]=0;
      vz[i]=0;
    }
  
  
  /* 3. start loop over time steps */ 
  
  for (int t=0; t<TMAX; t++)
    { 
      /* 4. call cicInterpolate() */ 
      
      cicInterpolate(ngrid, npart, x, y, z, rho);
      
      /* 5. call solvePoisson() */ 
      
      solvePoisson(a, rho, phi, ngrid);
      
      /* 6. write out data*/ 
      // to file?
      
      /* 7. call updateParticles() */ 
      updateParticles(ngrid, npart, a, da,&x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], &phi[0][0][0]);
      
      /* end loop */ 
    }
}
/*
 *Calculate contribution of each particle to the density grid points using the CIC interpolation.
 *Each particle contributes to the cell it is in and the neighboring grid cells based on the
 *algorithm in Section 2.8 of the write-up.
 */ 

void cicInterpolate(int ngrid, int npart, double *x, double *y, double *z, double ***rho) { 
    //declare parent cell locations
    int pcx, pcy, pcz;
    
    //declare distances from particle to cell (note that t_x/y/z is not explicitly defined but only used as '1 - d')
    int dx, dy, dz;
    
    //declare the factors that will be used for mass assignment (either dx or tx)
    int xfactor, yfactor, zfactor;
    
    //loop over particles
    for (int counter=0; counter<npart; counter++) {
        //get parent cell locations
        pcx = floor(x[counter]);
        pcy = floor(y[counter]);
        pcz = floor(z[counter]);
        
        //loop over relevant cells
        for (int i=pcx; i<pcx+2; i++) {
            //get the d and t variable for x
            dx = x[counter] - i;
            //if on the first iteration of x set the xfactor to tx, on the second iteration to dx
            if ((pcx - i) % 2 == 0) {xfactor = 1 - dx;}
            else {xfactor = dx;}
            for (int j=pcy; j<pcy+2; j++) {
                //get the d and t variable for y
                dy = y[counter] - j;
                //if on the first iteration of y set the yfactor to ty, on the second iteration to dy
                if ((pcy - i) % 2 == 0) {yfactor = 1 - dy;}
                else {yfactor = dy;}
                for (int k=pcz; k<pcz+2; k++) {
                    //get the d and t variable for z
                    dz = z[counter] - k;
                    //if on the first iteration of z set the zfactor to tz, on the second iteration to dz
                    if ((pcz - i) % 2 == 0) {zfactor = 1 - dz;}
                    else {zfactor = dz;}
                    
                    //now increment the mass of the appropriate rho index (assuming mass = 1)
                    rho[i % ngrid][j % ngrid][k % ngrid] += xfactor*yfactor*zfactor;
                }
            }
        }
        //now switch to the next particle
        counter++;
    }
}

/*  Solve poisson's equations and calculate acceleration field  */ 
void solvePoisson(double a, double ***rho, double ***phi, int ngrid) {

  double G[ngrid][ngrid][ngrid];      // Green's function
  const double G_Const= -( 3*OMEGA/(8*a) ); // constant in front of greens

  // start with rho(i, j, k), loop over all space
  // FFT rho to fourier space

  rho_fft = fftw3(rho); // loops internally

  // calulcate G
    // loop over cells l,m,n (In Fourier space!)
  for ( int l=-(0.5*ngrid); l<(0.5*ngrid + 1); l++ )
  {
    for ( int m=-(0.5*ngrid); m<(0.5*ngrid + 1); m++ )
    {
      for ( int n=-(0.5*ngrid); n<(0.5*ngrid + 1); n++ )
      {
        G[l][m][n] = pow(sin(pi*l/L), 2) + pow(sin(pi*m/L), 2) + pow(sin(pi*n/L), 2);
        G[l][m][n] = G_Const * pow(G[l][m][n], -1);

        phi[l][m][n] = G[l][m][n] * rho_fft[l][m][n];
      }
    }
  }

  // transform back into real space
  phi = fftw3(phi);


}

//Some helper functions for calculating the g for the x, y, and z directions 
double getGx(int i, int j, int k, vector<double> &phi){
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
}

double getGy(int i, int j, int k, vector<double> &phi){
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
}

double getGz(int i, int j, int k, vector<double> &phi){
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
}


/* Update position, velocities for each particle */ 
void updateParticles(int ngrid, int npart, double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, vector<double> &phi)
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
    gx = getGx(i, j, k, phi) * tx * ty * tz + getGx(i+1, j, k, phi) * dx * ty * tz + getGx(i, j+1, k, phi) * tx * dy * tz + getGx(i + 1, j+1, k. phi) * dx * dy * tz + getGx(i, j, k+1, phi)* tx * ty * dz + getGx(i+1, j, k+1, phi) * dx * ty * dz + getGx(i, j+1, k+1, phi) * tx * dy * dz + getGx(i + 1, j+1, k+1, phi) * dx * dy * dz;
    gy = getGy(i, j, k, phi) * tx * ty * tz + getGy(i+1, j, k, phi) * dx * ty * tz + getGy(i, j+1, k, phi) * tx * dy * tz + getGy(i + 1, j+1, k. phi) * dx * dy * tz + getGy(i, j, k+1, phi)* tx * ty * dz + getGy(i+1, j, k+1, phi) * dx * ty * dz + getGy(i, j+1, k+1, phi) * tx * dy * dz + getGy(i + 1, j+1, k+1, phi) * dx * dy * dz;
    gz = getGz(i, j, k, phi) * tx * ty * tz + getGz(i+1, j, k, phi) * dx * ty * tz + getGz(i, j+1, k, phi) * tx * dy * tz + getGz(i + 1, j+1, k. phi) * dx * dy * tz + getGz(i, j, k+1, phi)* tx * ty * dz + getGz(i+1, j, k+1, phi) * dx * ty * dz + getGz(i, j+1, k+1, phi) * tx * dy * dz + getGz(i + 1, j+1, k+1, phi) * dx * dy * dz;

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

    a = a + da; //Updates a accordingly
  }
}

int f(double a){
  //Needs to be implemented as shown onbottom of p8
  return pow((1 / a) * sqrt( OmegaM + OmegaK * a + OmegaL * pow(a, 2) ), -0.5);
}
