#include "cosmology.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    cout << setprecision(3) << fixed;
    
    // timing
    clock_t time_req = clock();
    
    // input parameters
    cosmology C;
    if (argc > 1) {
        C.OmegaM = atof(argv[1]);
        C.OmegaB = atof(argv[2]);
        C.zeq = atof(argv[3]);
        C.sigma8 = atof(argv[4]);
        C.h = atof(argv[5]);
        C.T0 = atof(argv[6]);
        C.ns = atof(argv[7]);
    }
    // run with PDG values
    else {
        C.OmegaM = 0.315;
        C.OmegaB = 0.0493;
        C.zeq = 3402.0;
        C.sigma8 = 0.811;
        C.h = 0.674;
        C.T0 = 2.7255;
        C.ns = 0.965;
    }
    
    C.initialize();
    cout << "The age of the universe today is " << C.age(0.0) << " Myr." << endl;
    cout << "The luminosity distance of a source at z = 1 is " << C.DL(1.0) << " kpc." << endl;
    
    // accuracy parameters
    int Nk = 1000;
    int NM = 1000;
    int Nz = 1000;
    
    // mass and redshift ranges
    double Mmin = 1.0, Mmax = 1.0e16;
    double zmin = 0.01, zmax = 20.0;

    // compute the variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    vector<vector<double> > sigma = C.sigmalist(Nk, NM, Mmin, Mmax);
    
    // compute the Seth-Tormen halo mass function {z,M,dndlnM}
    vector<vector<vector<double> > > dndlnM = C.hmflist(sigma, Nz, zmin, zmax);
    
    // output the halo mass function
    string filename = "hmf.dat";;
    ofstream outfile;
    outfile.open(filename.c_str());
    double z, M, hmf;
    for (int jz = 0; jz < Nz; jz++) {
        for (int jM = 0; jM < NM; jM++) {
            z = dndlnM[jz][jM][0];
            M = dndlnM[jz][jM][1];
            hmf = dndlnM[jz][jM][2];
            if (hmf < 1e-64) {
                hmf = 0.0;
            }
            outfile << z << "   " << M << "   " << hmf << endl;
        }
    }
    outfile.close();
    
    // random number generator
    rgen mt(time(NULL));
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
