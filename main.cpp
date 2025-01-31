#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"


// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    cout << setprecision(15) << fixed;
    
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
    vector<double> Zlist {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    
    // lensing source radius in kpc
    double rS = 10.0;
    
    // accuracy parameters for lensing
    int Nkappa = 40;
    int Nreal = 20000000;
    int Nbins = 2000;
    
    // read or generate lensing amplification distribution
    vector<vector<vector<double> > > Plnmuz = getPlnmu(C, Zlist, rS, Nkappa, Nreal, Nbins);
        
    // read UV luminosity data files, {MUV, Phi, +sigmaPhi, -sigmaPhi}
    vector<string> datafiles {"UVLF_ 2102.07775.txt", "UVLF_ 2108.01090.txt"};
    vector<vector<double> > data = readUVdata(datafiles);
    cout << "Read " << data.size() << " UV luminosity data points." << endl;
    
    cout << "Computing UV luminosity fit..." << endl;
    
    vector<double> initial {2.0e11, 0.075, -1.0, 0.30}; // first element in the MCMC chain, {M_c, epsilon, alpha, beta}
    vector<double> steps {0.2e11, 0.02, 0.04, 0.04}; // random walk step sizes
    vector<vector<double> > priors {{1.0e10, 1.0e12}, {0.02, 1.0}, {-1.6, -0.4}, {0.02, 1.2}}; // priors
    int Nsteps = 200; // number of steps
    int Nchains = 200; // number of chains
    rgen mt(time(NULL)); // random number generator
    
    vector<vector<double> > chain;
    vector<vector<vector<double> > > chains(Nchains);
    
    // output the MCMC chains and find the best fit
    outfile.open("MCMCchain.dat");
    int jmax;
    double logLmax = 0.0;
    vector<double> bf;
    for (int j = 0; j < Nchains; j++) {
        cout << "Generating chain number "<< j << "..." << endl;
        chain = mcmc_sampling(C, Zlist, Plnmuz, data, initial, steps, priors, Nsteps, mt);
        
        // find best fit
        jmax = max_element(chain.begin(), chain.end(), [](const vector<double> &a, const vector<double> &b) { return a.back() < b.back(); }) - chain.begin();
        if (chain[jmax][4] > logLmax) {
            logLmax = chain[jmax][4];
            bf = chain[jmax];
        }

        for (int j = 0; j < chain.size(); j++) {
            outfile << chain[j][0] << "   " << chain[j][1] << "   " << chain[j][2] << "   " << chain[j][3] << "   " << chain[j][4] << endl;
        }
    }
    outfile.close();
    
    cout << bf[0] << "   " << bf[1] << "   " << bf[2] << "   " << bf[3] << "   " <<bf[4] << endl;
    
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
