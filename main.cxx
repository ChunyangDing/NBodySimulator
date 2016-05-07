#include <fftw3.h>
#include <math.h> // for floor(), pow()
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h> // for malloc, drand48()

using namespace std;

// function prototypes
void cicInterpolate(int ngrid, int npart, double *x, double *y, double *z, vector<double>& rho);
void solvePoisson(double a, vector<double>& rho,  vector<double>& phi, int ngrid);
void updateParticles(int ngrid, int npart, double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, double ***phi);
int f(double a);

void usage(const char* prog) {
    cerr << "Usage: " << prog << " [-V] [-h] [fileName]" << endl;
}

int main(int argc, char* argv[]) {
  int ch;
  //update code version if any changes are made to the code
  char version[] = "v.74";
  
  while ((ch = getopt(argc, argv,"Vh")) != -1) {
        //note that there break statements are redundant for cases when the program is quitted
    switch (ch) {
      case 'V':
        cout << "Version of code: " << version << endl;
        //quit the program if 'V' is passed.
        return 0;
      case 'h':
        cout << "This program will simulate dark matter based on given input conditions and constants.  It will create a text file at the end with the coordinates and velocities of each particle at each time-step.  Please call the program in the following way:" << endl;
        usage(argv[0]);
        cout << "All the arguments in the program are optional.\n-V\tprints the program version and quits\n-h\tprints instructions for the program and quits.\n[fileName]\t is an optional argument to specify the name of the output file." << endl;
        //quit the program if 'h' is passed.
        return 0;
      case '?':
        //call usage function if unkown option
        cout << "You have passed an unspecified option.  Run code as '" << argv[0] << " -h' for help." << endl;
				usage(argv[0]);
        //return 2 and quit to specify that the program has did not recieve proper options
        return 2;
      }
  }
    
  //now check for the file name
  int optLeftOver = argc - optind;
  if (optLeftOver > 1) {
    cout << "You have entered more than one standard arguments or did not put [file] at the end.  Run code as '" << argv[0] << " -h' for help." << endl;
    //return 2 and quit to specify that the program has did not recieve proper argument
    return 2;
  }
    
  //declare the output file
  ofstream outFile;
  //if there is one opt left then that will be the name, else the file will be named by default
  if (optLeftOver == 1) {
    outFile.open(argv[optind]);
  }
  else {outFile.open("out.txt");}


  // relevant variables
  const int ngrid = pow(64,1);
  const int npart = pow(32,1);

  // cosmological variables
  const double OmegaM = 0.27;
  const double OmegaL = 0.73;
  const double OmegaK = 0.0000824;
  const double G = 6.6740831 * pow(10, -11);
  const double rhoCrit = 1.06 * pow(10, -29);
  const double H = 67.80;

  const double pi = 3.141592654;

  // variables for converting to actual units (we can ignore these for our purposes)
  //const double r0 = 1;
  //const double t0 = 1/H;
  //const double v0 = r0 / t0;
  //const double rho0 = ((3 * pow(H, 2)) / (8 * pi * G)) * OmegaM;
  //const double phi0 = pow(v0, 2);
  const double L = 100.0; // scales the size of the model

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
  
  vector<double> rho;    // Should describe mass density for each cell
  vector<double> phi;     // Should have unique density for each cell
  

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

  
    
  // setup initial conditions
  for (int i=0; i<ngrid; i++)
  {
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
    x[i] = x[i] + sin( 2*pi*x[i] / L);
    y[i] = y[i] + sin( 2*pi*y[i] / L);
    z[i] = z[i] + sin( 2*pi*z[i] / L);

    vx[i]=0;  // initial velocities are all zero
    vy[i]=0;  
    vz[i]=0;  
  }
  
  //Make a header for the outputfile
  outFile << "ID\tx\ty\tz\tvx\tvy\tvz" << endl;
  
  // loop over time-steps
  while ( a < aMAX )
    { 
      /* 4. call cicInterpolate() */ 
      
      cicInterpolate(ngrid, npart, x, y, z, rho);
      
      /* 5. call solvePoisson() */ 
      
      solvePoisson(a, rho, phi, ngrid);
      
      /* 6. write out data*/
      outFile << "\nAt a = " << a;
      for (int i=0; i<npart; i++) {
        outFile << "\n" << x[i] << " "  << y[i] << " "  << z[i] << " "  << vx[i] << " "  << vy[i] << " "  << vz[i];
      }
      
      /* 7. call updateParticles() */ 
      updateParticles(ngrid, npart, a, da,&x[0], &y[0], &z[0], &vx[0], &vy[0], &vz[0], phi);
      
      /* end loop */ 
      a = a + da;
    }

    return 0;
}
/*
 *Calculate contribution of each particle to the density grid points using the CIC interpolation.
 *Each particle contributes to the cell it is in and the neighboring grid cells based on the
 *algorithm in Section 2.8 of the write-up.
 */ 
void cicInterpolate(int ngrid, int npart, double *x, double *y, double *z, vector<double>& rho) {
    //declare parent cell locations
    int pcx, pcy, pcz;
    
    //declare distances from particle to cell (note that t_x/y/z is not explicitly defined but only used as '1 - d_x/y/x')
    double dx, dy, dz;
    
    //declare the factors that will be used for mass assignment (either d or t)
    double xfactor, yfactor, zfactor;
    
    //loop over particles
    for (int counter=0; counter<npart; counter++) {
        //get parent cell locations
        pcx = floor(x[counter]);
        pcy = floor(y[counter]);
        pcz = floor(z[counter]);
        
        dx = x[counter] - pcx;
        dy = y[counter] - pcy;
        dz = z[counter] - pcz;
        
        //loop over relevant cells
        for (int i=pcx; i<pcx+2; i++) {
            //get the d and t variable for x
            //if on the first iteration of x set the xfactor to tx, on the second iteration to dx
            if (i - pcx == 0) {xfactor = 1 - dx;}
            else {xfactor = dx;}
            for (int j=pcy; j<pcy+2; j++) {
                //get the d and t variable for y
                //if on the first iteration of y set the yfactor to ty, on the second iteration to dy
                if (j - pcy == 0) {yfactor = 1 - dy;}
                else {yfactor = dy;}
                for (int k=pcz; k<pcz+2; k++) {
                    //get the d and t variable for z
                    //if on the first iteration of z set the zfactor to tz, on the second iteration to dz
                    if (k - pcz == 0) {zfactor = 1 - dz;}
                    else {zfactor = dz;}
                    
                    //now increment the mass of the appropriate rho index (assuming mass = 1)
                    rho[(i%ngrid)*ngrid*ngrid + (j%ngrid)*ngrid + (k%ngrid)] += xfactor*yfactor*zfactor;
                }
            }
        }
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
