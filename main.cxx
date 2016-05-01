/* Include any necessary .h files, such as math.h, stdio.h, etc */ 
#include <fftw3.h>
#include <math.h> // for floor(), pow()

using namespace std;

int main() 
{
  /* 1. declarations of relevant variables */
  const int ngrid=pow(64,1); // supposed to be 64**3, but don't want to break things with big numbers...
  const int npart=pow(32,1); // supposed to be 32**3, but big numbers are bad when you make mistakes...
  
  const int TMAX=100; // number of timesteps

  const double OmegaM = 0.27;
  const double OmegaL = 0.73;
  const double OmegaK = 0.0000824;
  const double G = 6.6740831 * pow(10, -11);
  const double rhoCrit = 1.06 * pow(10, -29);
  const double H = 67.80;
  const double pi = 3.141592654;

  const double r0 = 1;
  const double t0 = 1/H;
  const double v0 = r0 / t0;
  const double rho0 = ((3 * pow(H, 2)) / (8 * M_PI * G)) * OmegaM;
  const double phi0 = pow(v0, 2);

  double x[npart]; //Should describe every particle
  double y[npart];
  double z[npart];
  
  double vx[npart];
  double vy[npart];
  double vz[npart];
  
  double rho[ngrid][ngrid][ngrid];    // Should describe mass density for each cell??
  double phi[ngrid][ngrid][ngrid];     // Should have unique density for each cell
  
  double a;     // Should be the scale parameter, usually = 1
  double da;    // Some small change correlated with time step.
  
  
  //There are a lot of variables that are used throughout here...
  
  /* 2. read in or setup initial conditions */
  
  // I think this part is explained in the hand-out
  // I have it set up in loops write now just so there is something in the arrays
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

/* Update position, velocities for each particle */ 
void updateParticles(int ngrid, int npart, double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, double ***phi)
{
  for (int abc = 0; abc < npart; abc++){
    /* 1. declarations of relevant variables */ 
    int i = static_cast<int>(floor(*x));
    int j = static_cast<int>(floor(*y));
    int k = static_cast<int>(floor(*z));
    double gx = 0;
    double gy = 0;
    double gz = 0;

    /* 2. calculate particle accelerations from phi */
    if (i == 0){
      gx = -(*(*(*(phi + (i)) + j ) + k) - *(*(*(phi + (i + 1))+j)+k))/2.0;
    }
    else {
      if (i == ngrid){
	gx = -(*(*(*(phi + (i - 1)) + j ) + k) - *(*(*(phi + (i))+j)+k))/2.0;      
      }
      else{
	gx = -(*(*(*(phi + (i - 1)) + j ) + k) - *(*(*(phi + (i + 1))+j)+k))/2.0;
      }
    }

    if (j == 0){
      gy = -(*(*(*(phi + i) + (j+1)) + k ) - *(*(*(phi + i) +(j)) + k) ) /2.0;
    }
    else {
      if (j == ngrid){
	gy = -(*(*(*(phi + i) + (j)) + k ) - *(*(*(phi + i) +(j-1)) + k) ) /2.0;
      }
      else{
	gy = -(*(*(*(phi + i) + (j+1)) + k ) - *(*(*(phi + i) +(j-1)) + k) ) /2.0;
      }
    }

    if (k == 0){
      gz = -(*(*(*(phi + i) + j) + (k+1)) - *(*(*(phi + i) + j) + (k) ) )/2.0;
    }
    else {
      if (k == ngrid){
	gz = -(*(*(*(phi + i) + j) + (k)) - *(*(*(phi + i) + j) + (k-1) ) )/2.0;
      }
      else {
	gz = -(*(*(*(phi + i) + j) + (k+1)) - *(*(*(phi + i) + j) + (k-1) ) )/2.0;
      }
    }
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

int f(double a){
  //Needs to be implemented as shown onbottom of p8
  return pow((1 / a) * sqrt( OmegaM + OmegaK * a + OmegaL * pow(a, 2) ), -0.5)
}
