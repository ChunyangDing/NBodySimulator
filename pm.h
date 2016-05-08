#include <vector>

void cicInterpolate(double *x, double *y, double *z, std::vector<double>& rho);
void solvePoisson(double a, std::vector<double> &rho, fftw_complex *frho,fftw_complex *fphi, std::vector<double> &phi);
double getGx(int i, int j, int k, std::vector<double> &phi);
double getGy(int i, int j, int k, std::vector<double> &phi);
double getGz(int i, int j, int k, std::vector<double> &phi);
void updateParticles(double a, double da, double *x, double *y, double *z, double *vx, double *vy, double *vz, std::vector<double> &phi);
double f(double a);
void usage(const char* prog);

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

const double L = 100.0; // scales the size of the model

// variables for converting to actual units (we can ignore these for our purposes)
//const double r0 = 1;
//const double t0 = 1/H;
//const double v0 = r0 / t0;
//const double rho0 = ((3 * pow(H, 2)) / (8 * pi * G)) * OmegaM;
//const double phi0 = pow(v0, 2);