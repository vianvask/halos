
/*
 ./uvlf 1 0 4 15.0 &
 ./uvlf 1 0 4 30.0 &
 ./uvlf 1 0 5 15.0 &
 ./uvlf 1 0 5 30.0
 
 */

#include "cosmology.h"
#include "UVluminosity.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(6) << fixed;
    
    const int doUVfit = atoi(argv[1]); // 0: no, 1: yes
    const int pdm = atoi(argv[2]); // 0: linear priors, 1: log priors for the beyond LCDM parameter
    const int dm = atoi(argv[3]); // 0: cold DM, 1: fuzzy DM, 2: warm DM, 3: white noise, 4: infl. magnetic fields, 5: PT magnetic fields
    const double zcut = atof(argv[4]); // redshift cut on uvlf data: 15.0 removes the uncertain high-z data, 30.0 doesn't remove anything

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
    C.NM = 240;
    
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
    
    C.outdir = "dataU";
    C.initialize(dm);
    
    cout << "Reading lensing amplifications..." << endl;
    C.Zlist = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5, 17.0, 25.0};
    C.Plnmuz = getPlnmu(C.outdir/"Plnmuz.dat");
    if (C.Plnmuz.size() != C.Zlist.size()) {
        cout << "Wrong Plnmu file." << endl;
        return 0;
    }
    
    // remove z > z_cut
    C.Zlist.erase(remove_if(C.Zlist.begin(), C.Zlist.end(), [zcut](double x) { return x > zcut; }), C.Zlist.end() );
    
    cout << "Computing UV luminosity functions..." << endl;
   
    // parameters: (logMt, Mc, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0, sigmaUV, logm)
    vector<double> bf = {8.0, 3.78e11, 0.062, 1.0, 0.4, 0.1, 10.5, 0.37, 11.5, 23.0, 0.08, -1.3};
    
    if (dm == 4) {
        bf = {8.0, 3.78e11, 0.062, 1.0, 0.4, 0.1, 10.5, 0.37, 11.5, 23.0, 0.12, 0.5};
    }
    
    if (dm == 5) {
        bf = {8.0, 3.78e11, 0.062, 1.0, 0.4, 0.1, 10.5, 0.37, 11.5, 23.0, 0.12, 0.003};
    }
    
    if (doUVfit == 1) {
        int Nsteps = 40000; // chain length without burn-in
        int Nburnin = 2000; // burn-in
        int Nchains = 4; // number of chains
        double xstep = 20.0; // step size = prior range/xstep
        
        // parameters: (logMt, Mc, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0, sigmaUV)
        vector<vector<double> > priors = {{6.5 ,9.5}, {2.4e11, 5.2e11}, {0.0562, 0.0690}, {0.65, 1.25}, {0.2, 0.6}, {0.05, 0.6}, {9.8, 11.4}, {0.05, 0.7}, {7.0, 25.0}, {16.0, 36.0}, {0.08, 0.24}};
        
        if (dm == 1) {
            if (pdm == 0) {
                priors.push_back({C.m22list.front(), C.m22list.back()});
            } else {
                priors.push_back({log10(C.m22list.front()), log10(C.m22list.back())});
            }
        }
        if (dm == 2) {
            if (pdm == 0) {
                priors.push_back({C.m3list.front(), C.m3list.back()});
            } else {
                priors.push_back({log10(C.m3list.front()), log10(C.m3list.back())});
            }
        }
        if (dm == 3) {
            if (pdm == 0) {
                priors.push_back({C.kclist.front(), C.kclist.back()});
            } else {
                priors.push_back({log10(C.kclist.front()), log10(C.kclist.back())});
            }
        }
        if (dm == 4 || dm == 5) {
            if (pdm == 0) {
                priors.push_back({C.Blist.front(), C.Blist.back()});
            } else {
                priors.push_back({log10(C.Blist.front()), log10(C.Blist.back())});
            }
            
            // extend prior of alpha, beta and sigma
            priors[3] = {0.65, 1.45};
            priors[10] = {0.1, 0.24};
        }
        
        bf = UVLFfit(C, priors, Nsteps, Nburnin, Nchains, xstep, dm, pdm);
    }
    
    // output the UV luminosity function for the best fit
    writeUVLF(C, bf[0], bf[1], bf[2], bf[3], bf[4], bf[5], bf[6], bf[7], bf[8], bf[9], bf[10], bf[11], dm, pdm);
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
