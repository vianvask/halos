#include "cosmology.h"
#include "lensing.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(8) << fixed;
    
    // cosmological parameters: PDG values
    cosmology C;
    C.OmegaM = 0.315;
    C.OmegaB = 0.0493;
    C.zeq = 3402.0;
    C.sigma8 = 0.811;
    C.h = 0.674;
    C.T0 = 2.7255;
    C.ns = 0.965;
    
    // halo masses
    C.Mmin = 1.0e7;
    C.Mmax = 1.0e17;
    C.NM = 100;
    
    // redshifts
    C.zmin = 0.01;
    C.zmax = 16.01;
    C.Nz = 160;
    
    cout << "Computing halo mass functions..." << endl;
    C.initialize(0);
    
    cout << "Generating lensing amplifications..." << endl;
    
    int Nreal = 1e6; // realizations
    int Nbins = 200; // P(lnmu) bins
    double rS = 0.0; // lensing source radius in kpc
    
    rgen mt(time(NULL)); // random number generator
    
    // list of redshift values at which P(lnmu) is computed
    //C.Zlist = {1.0};
    C.Zlist = {0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
    //C.Zlist = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5, 17.0, 25.0};
    C.Plnmuz = getPlnmu(C, rS, Nreal, Nbins);
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
