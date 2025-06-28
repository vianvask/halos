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
    
    rgen mt(time(NULL)); // random number generator
    
    cout << "Generating lensing amplifications..." << endl;
    
    int Nhalos = 10;
    int Nreal = 2e4; // realizations
    int Nbins = 12; // P(lnmu) bins
    double rS = 0.0; // lensing source radius in kpc
    
    // list of redshift values at which P(lnmu) is computed
    //vector<double> Zlist = {10.0};
    vector<double> Zlist = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
    // vector<double> Zlist = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5, 17.0, 25.0};
    
    ofstream outfile;
    outfile.open("Plnmu.dat");
    vector<vector<double> > Plnmu;
    for (int jz = 0; jz < Zlist.size(); jz++) {
        cout << "z = " << Zlist[jz] << endl;
        Plnmu = Plnmuf(C, Zlist[jz], rS, Nhalos, Nreal, Nbins, mt);
        for (int jb = 0; jb < Plnmu.size(); jb++) {
            outfile << Zlist[jz] << "   " << Plnmu[jb][0] << "   " << Plnmu[jb][1] << endl;
        }
    }
    outfile.close();
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
