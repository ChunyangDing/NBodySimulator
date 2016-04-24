/* Include any necessary .h files, such as math.h, stdio.h, etc */ #include <fftw3.h> 
int main() 
{
  /* 1. declarations of relevant variables */
  /* 2. read in or setup initial conditions */
  /* 3. start loop over time steps */ 
  /* 4. call cicInterpolate() */ 
  /* 5. call solvePoisson() */ 
  /* 6. write out data*/ 
  /* 7. call updateParticles() */ 
  /* end loop */ 
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
