void solvePoisson(double a, double *rho, double *phi) {
    /* 1. declarations of relevant variables */ 
    int N=ngrid;
	fftw_complex *frho= fftw_malloc(sizeof(fftw_complex)*N*N*(N/2 +1));
	fftw_plan pf;
	fftw_plan pb;
	double kx, ky, kz;

	fftw_complex *fphi= fftw_malloc(sizeof(fftw_complex)*N*N*(N/2 + 1));

	double *green= new double[N*N*(N/2 +1)];

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
    	xcounter++;
    for (int m=-(N/2); m<(N/2); m++){
    	ycounter++;        	
    for (int n=0; n<(N/2 + 1); n++){
    	zcounter++;
    	if ( (l==0) && (m==0) && (n==0)) green[xcounter*N*(N/2 +1) + ycounter*(N/2 +1) + zcounter] = 0;
    	kx = 2.0*pi*l/N;
    	ky = 2.0*pi*m/N;
    	kz = 2.0*pi*n/N;
    	// determine green
    	green[xcounter*N*(N/2 +1) + ycounter*(N/2 +1) + zcounter] = -((3.0*OmegaM)/(8.0*a))*pow( pow( sin(.5*kx),2) + sin(.5*ky),2) + sin(.5*kz),2),-1);
    }}}

    for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
    for (int k=0; k<(N/2 +1); k++){
    	fphi[i*N*(N/2 +1) + j*(N/2 + 1) + k] = green[i*N*(N/2 +1) + j*(N/2 + 1) + k] * frho[i*N*(N/2 +1) + j*(N/2 + 1) + k];
    }}}


    /* 4. reverse transformation using fftw_plan_dft_c2r_3d() and fftw_execute() */
    pb = fftw_plan_dft_c2r_3d(ngrid, ngrid, ngrid, fphi, phi, FFTW_ESTIMATE);

    fftw_execute(pb);

    for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
    for (int k=0; k<N; k++){
    	phi[i*N*N + j*N + k] *= 1.0/(N*N*N);
    }}} 

    delete[] green;
    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb);
    fftw_free(frho);
    fftw_free(fphi);
} 