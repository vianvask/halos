#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"


// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(14) << fixed;
    
    cout << "Preparing..." << endl;
    cosmology C;
    
    // cosmological parameters: PDG values
    C.OmegaM = 0.315;
    C.OmegaB = 0.0493;
    C.zeq = 3402.0;
    C.sigma8 = 0.811;
    C.h = 0.674;
    C.T0 = 2.7255;
    C.ns = 0.965;
    
    // mass and redshift ranges
    C.Mmin = 1.0e6;
    C.Mmax = 1.0e17;
    C.zmin = 0.01;
    C.zmax = 20.01;
    
    C.Nk = 1000; // integral over k in computation of sigma
    C.NM = 1000; // number of halo mass bins
    C.Nz = 100; // number of redshift bins
    
    double z, M, HMF, dotM, DdotM;
    ofstream outfile;
    
    // output some halo mass functions
    vector<double> m22list {0.0, 1.0, 10.0, 100.0};
    for (double m22 : m22list) {
        C.initialize(m22);
        outfile.open("HMF_m22_" + to_string_prec(m22,1) + ".dat");
        for (int jz = 0; jz < C.Nz; jz++) {
            for (int jM = 0; jM < C.NM; jM++) {
                z = C.zlist[jz];
                M = C.HMFlist[jz][jM][1];
                HMF = C.HMFlist[jz][jM][2];
                dotM = C.dotMlist[jz][jM][2];
                DdotM = C.dotMlist[jz][jM][3];
                
                outfile << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << "   " << DdotM << endl;
            }
        }
        outfile.close();
    }

    C.initialize(1.0);
    
    // read or generate lensing amplification distribution
    int Nkappa = 40;
    int Nreal = 10000000;
    int Nbins = 200;
    double rS = 10.0;  // lensing source radius in kpc
    vector<double> Zlist {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5};
    vector<vector<vector<double> > > Plnmuz = getPlnmu(C, Zlist, rS, Nkappa, Nreal, Nbins);
    
    // read UV luminosity data files, {MUV, Phi, +sigmaPhi, -sigmaPhi}
    vector<string> datafiles {"UVLF_2102.07775.txt", "UVLF_2108.01090.txt", "UVLF_2403.03171.txt"};
    vector<vector<double> > data = readUVdata(datafiles);
    cout << "Read " << data.size() << " UV luminosity data points." << endl;
    
    cout << "Computing UV luminosity fit..." << endl;
    
    // MCMC over parameters {M_c, epsilon, alpha, beta}
    int Nsteps = 100; // number of steps
    int Nchains = 40; // number of chains
    vector<vector<double> > priors {{2.0e11, 7.0e11}, {0.02, 0.07}, {0.7, 1.5}, {0.1, 0.7}, {0.5, 2.0}, {7.0, 11.0}}; // flat priors
    vector<double> initial(6,0.0); // first element
    vector<double> steps {0.2e11, 0.02, 0.05, 0.05, 0.05, 0.1}; // random walk step sizes
    rgen mt(time(NULL)); // random number generator

    // output the MCMC chains and find the best fit
    outfile.open("MCMCchains.dat");
    int jmax;
    double logLmax = 0.0;
    vector<double> bf;
    vector<vector<double> > chain;
    for (int j = 0; j < Nchains; j++) {
        cout << j << endl;
        
        // generate initial point inside the prior ranges
        for (int jp = 0; jp <  initial.size(); jp++) {
            initial[jp] = randomreal(priors[jp][0],priors[jp][1]);
        }
        
        // generate MCMC chain
        chain = mcmc_sampling(C, Zlist, Plnmuz, data, initial, steps, priors, Nsteps, mt);
        
        // find best fit
        jmax = max_element(chain.begin(), chain.end(), [](const vector<double> &a, const vector<double> &b) { return a.back() < b.back(); }) - chain.begin();
        if (chain[jmax].back() > logLmax) {
            logLmax = chain[jmax].back();
            bf = chain[jmax];
        }

        // output the MCMC chain
        for (int j = 0; j < chain.size(); j++) {
            for (int jp = 0; jp < chain[0].size(); jp++) {
                outfile << chain[j][jp] << "   ";
            }
            outfile << endl;
        }
    }
    outfile.close();
    
    // output the UV luminosity function for the best fit
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, Zlist, Plnmuz, 1.0e9, bf[0], bf[1], bf[2], bf[3], bf[4], bf[5]);
    outfile.open("UVluminosity.dat");
    for (int jZ = 0; jZ < Zlist.size(); jZ++) {
        z = Zlist[jZ];
        for (int jM = 0; jM < C.NM; jM++) {
            // {z, M, MUV, no dust + no lensing, dust + no lensing, dust + lensing}
            outfile << z << "   " << PhiUVlist[jZ][jM][0] << "   " << PhiUVlist[jZ][jM][4] << "   " << max(1.0e-64,PhiUVlist[jZ][jM][1]) << "   " << max(1.0e-64,PhiUVlist[jZ][jM][2]) << "   " << max(1.0e-64,PhiUVlist[jZ][jM][3]) << endl;
        }
    }
    outfile.close();

    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
