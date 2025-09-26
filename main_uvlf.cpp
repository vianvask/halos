#include "cosmology.h"
#include "UVluminosity.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(6) << fixed;
    
    const int doUVfit = atoi(argv[1]); // 0: no, 1: yes
    const int dm = atoi(argv[2]); // 0: cold DM, 1: fuzzy DM, 2: warm DM, 3: white noise
        
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
    C.Mmin = 1.0e6;
    C.Mmax = 1.0e17;
    C.NM = 2000;
    
    // redshifts
    C.zmin = 0.01;
    C.zmax = 36.01;
    C.Nz = 180;

    // FDM masses in 10^-22 eV
    C.m22min = 1.0;
    C.m22max = 2000.0;
    C.Nm22 = 50;
    
    // WDM masses in keV
    C.m3min = 0.7;
    C.m3max = 40.0;
    C.Nm3 = 50;
    
    // kc values in 1/kpc
    C.kcmin = 0.006;
    C.kcmax = 2.0;
    C.Nkc = 50;
    
    cout << "Computing halo mass functions..." << endl;
    C.initialize(dm);
    
    cout << "Reading lensing amplifications..." << endl;
    C.Zlist = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5, 17.0, 25.0};
    C.Plnmuz = getPlnmu();    
    if (C.Plnmuz.size() != C.Zlist.size()) {
        cout << "Wrong Plnmu.dat file." << endl;
        return 0;
    }
    
    cout << "Computing UV luminosity functions..." << endl;
    
    // parameters: (logMt, Mc, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0, sigmaUV, logm)
    vector<double> bf = {7.94, 4.00e11, 0.0612, 0.882, 0.404, 0.272, 10.73, 0.292, 12.7, 22.9, 0.068, -1.6};
    //vector<double> bf = {7.94, 4.00e11, 0.0612, 0.882, 0.404, 0.272, 10.73, 0.292, 12.7, 22.9, 0.068, 0.75};
    //vector<double> bf = {7.94, 4.00e11, 0.0612, 0.882, 0.404, 0.272, 10.73, 0.292, 26, 36, 0.068, 0.75};
    //vector<double> bf = {7.94, 4.00e11, 0.0612, 0.882, 0.404, 0.272, 26, 1.0, 26, 36, 0.068, 0.75};

    if (doUVfit == 1) {
        int Nsteps = 10000; // chain length without burn-in
        int Nburnin = 2000; // burn-in
        int Nchains = 8; // number of chains
        double xstep = 16.0; // step size = prior range/xstep
        
        // prior of the beyond CDM parameter is determined from the range given above
        vector<vector<double> > priors = {{6.0,10.0}, {3.0e11, 5.0e11}, {0.0563, 0.0658}, {0.6, 1.1}, {0.2, 0.6}, {0.05, 0.6}, {9.8, 11.4}, {0.05, 0.7}, {7.0, 25.0}, {16.0, 36.0}, {0.05, 0.2}};
        
        bf = UVLFfit(C, priors, Nsteps, Nburnin, Nchains, xstep, dm);
    }
    
    // output the UV luminosity function for the best fit
    writeUVLF(C, bf[0], bf[1], bf[2], bf[3], bf[4], bf[5], bf[6], bf[7], bf[8], bf[9], bf[10], bf[11], dm);

    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
