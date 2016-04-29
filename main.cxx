/* Include any necessary .h files, such as math.h, stdio.h, etc */ 
#include <fftw3.h>
#include <cmath> // for pow() 
int main() 
{
  /* 1. declarations of relevant variables */
  const int ngrid=pow(64,1); // supposed to be 64**3, but don't want to break things with big numbers...
  const int npart=pow(32,1); // supposed to be 32**3, but big numbers are bad when you make mistakes...

  const int TMAX=100; // number of timesteps

  double x[ngrid];
  double y[ngrid];
  double z[ngrid];

  double vx[ngrid];
  double vy[ngrid];
  double vz[ngrid];

  double rho[2][3][ngrid];    // I don't really know what is going on here.... at all.
  double phi[10][10][10];     // ditto

  double a;     // don't know what this is (scale parameter?)
  double da;    // time-stepping parameter?


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
    updateParticles(ngrid, npart, a, da, x, y, z, vx, vy, vz, phi);

    /* end loop */ 
  }
}
/*
 *Calculate contribution of each particle to the density grid points using the CIC interpolation.
 *Each particle contributes to the cell it is in and the neighboring grid cells based on the
 *algorithm in Section 2.8 of the write-up.
 */ 

void cicInterpolate(int ngrid, int npart, double *x, double *y, double *z, double ***rho) { 
  /* 1. declarations of relevant variables */ 
  /* 2. loop over particles and apply the CIC interpolation to calculate density of each cell */
}

/*  Solve poisson's equations and calculate acceleration field  */ 
void solvePoisson(double a, double ***rho, double ***phi, int ngrid) { 
  /* 1. declarations of relevant variables */ 
  /* 2. fourier transform rho into fourier space 
     We suggest looking to the fftw_plan_dft_r2c_3d() function for this transform. 
     Any fourier transform in fftw is then followed by the execution command.  
     For example: 
    // setup fftw plan 
    pf = fftw_plan_dft_r2c_3d(ngrid, ngrid, ngrid, in, out, FFTW_ESTIMATE); 
    // take fourier transform
    fftw_execute(pf); 
  */
  /* 3. calculate green's function in fourier space */ 
  /* 4. reverse transformation using fftw_plan_dft_c2r_3d() and fftw_execute() */
}

/* Update position, velocities for each particle */ 
void updateParticles(int ngrid, int npart, double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, double ***phi){     
  /* 1. declarations of relevant variables */ 
  /* 2. calculate particle accelerations from phi */ 
  /* 3. update particle velocities */ 
  /* 4. update particle positions */ 
}
