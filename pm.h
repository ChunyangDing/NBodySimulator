///////////////////////////////////////////////////////////////////////////////
//    pm.h :: header file
//
//	  not much to see here. Just global variables and function prototypes
//
///////////////////////////////////////////////////////////////////////////////


void cicInterpolate(double *x, double *y, double *z,double *rho);
void solvePoisson(double a, double *rho, double *phi);
void updateParticles(double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *phi);


void printVec3D(double* rho);


double f(double a);

const int N = 32*32;
const int M1D = 64;
const int Mtot = M1D*M1D*M1D;

const double OmegaM=0.27;
const double OmegaL = 0.73;
const double OmegaK = 0.0000824;

const double pi=3.14159;

double Dplus;
