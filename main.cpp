#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"


// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    cout << setprecision(4) << fixed;
    
    cout << "Preparing..." << endl;
    
    // timing
    clock_t time_req = clock();

    // input parameters: PDG values
    cosmology C;
    C.OmegaM = 0.315;
    C.OmegaB = 0.0493;
    C.zeq = 3402.0;
    C.sigma8 = 0.811;
    C.h = 0.674;
    C.T0 = 2.7255;
    C.ns = 0.965;
    
    // accuracy parameters
    C.Nk = 1000;
    C.NM = 1200;
    C.Nz = 100;
    
    // mass and redshift ranges
    C.Mmin = 1.0e5;
    C.Mmax = 1.0e17;
    C.zmin = 0.01;
    C.zmax = 20.01;
    
    C.initialize();
    
    double z, M, HMF, dotM, DdotM;
    ofstream outfile;
    
    // output the halo mass function and the UV luminosity function
    outfile.open("HMF.dat");
    for (int jz = 0; jz < C.Nz; jz++) {
        for (int jM = 0; jM < C.NM; jM++) {
            z = C.zlist[jz];
            M = C.HMFlist[jz][jM][1];
            HMF = C.HMFlist[jz][jM][2];
            dotM = C.dotMlist[jz][jM][2];
            DdotM = C.dotMlist[jz][jM][3];
            
            outfile << z << "   " << M << "   " << max(1.0e-64,HMF) << "   " << dotM << "   " << DdotM << endl;
        }
    }
    outfile.close();
    
    // list of z values for lensing
    vector<double> Zlist {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5};
    
    // lensing source radius in kpc
    double rS = 10.0;
    
    // accuracy parameters for lensing
    int Nkappa = 40;
    int Nreal = 10000000;
    int Nbins = 200;
    
    // read or generate lensing amplification distribution
    vector<vector<vector<double> > > Plnmuz = getPlnmu(C, Zlist, rS, Nkappa, Nreal, Nbins);
        
    // read UV luminosity data files, {MUV, Phi, +sigmaPhi, -sigmaPhi}
    vector<string> datafiles {"UVLF_2102.07775.txt", "UVLF_2108.01090.txt", "UVLF_2403.03171.txt"};
    vector<vector<double> > data = readUVdata(datafiles);
    cout << "Read " << data.size() << " UV luminosity data points." << endl;
    
    cout << "Computing UV luminosity fit..." << endl;
    
    // MCMC over parameters {M_c, epsilon, alpha, beta}
    int Nsteps = 200; // number of steps
    int Nchains = 200; // number of chains
    vector<vector<double> > priors {{0.6e11, 4.4e11}, {0.04, 0.11}, {-1.5, -0.5}, {0.1, 0.7}}; // flat priors
    vector<double> initial(4,0.0); // first element
    vector<double> steps {0.2e11, 0.02, 0.04, 0.04}; // random walk step sizes
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
        if (chain[jmax][4] > logLmax) {
            logLmax = chain[jmax][4];
            bf = chain[jmax];
        }

        // output the MCMC chain
        for (int j = 0; j < chain.size(); j++) {
            outfile << chain[j][0] << "   " << chain[j][1] << "   " << chain[j][2] << "   " << chain[j][3] << "   " << chain[j][4] << endl;
        }
    }
    outfile.close();
        
    // output the UV luminosity function (no dust + no lensing, dust + no lensing, dust + lensing) for the best fit
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, Zlist, Plnmuz, 1.0e9, bf[0], bf[1], bf[2], bf[3], 1.0);
    outfile.open("UVluminosity.dat");
    for (int jZ = 0; jZ < Zlist.size(); jZ++) {
        z = Zlist[jZ];
        for (int jM = 0; jM < C.NM; jM++) {
            outfile << z << "   " << PhiUVlist[jZ][jM][0] << "   " << max(1.0e-64,PhiUVlist[jZ][jM][1]) << "   " << max(1.0e-64,PhiUVlist[jZ][jM][2]) << "   " << max(1.0e-64,PhiUVlist[jZ][jM][3]) << endl;
        }
    }
    outfile.close();

    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
