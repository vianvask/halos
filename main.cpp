#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(4) << fixed;
    
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
    
    // integral over k in computation of sigma
    C.Nk = 1000;
    
    // halo masses
    C.Mmin = 1.0e5;
    C.Mmax = 1.0e17;
    C.NM = 1200;
    
    // redshifts
    C.zmin = 0.01;
    C.zmax = 20.01;
    C.Nz = 100;
    
    // FDM masses in 10^-22 eV
    C.m22min = 0.8;
    C.m22max = 80.0;
    C.Nm22 = 0;
    
    // WDM masses in keV
    C.m3min = 1.0;
    C.m3max = 100.0;
    C.Nm3 = 10;
    
    C.initialize();
    ofstream outfile;
    
    cout << "Computing halo mass functions (" << C.m22list.size() << " m22 values and " << C.m3list.size() << " m3 values)..." << endl;
    C.FDM_halos();
    C.WDM_halos();
    C.CDM_halos();
    
    // output the CDM halo mass function
    double m22, m3, z, M, sigma, dsigma, HMF, dotM;
    outfile.open("sigma_CDM.dat");
    for (int jM = 0; jM < C.NM; jM++) {
        M = C.sigmalist[jM][0];
        sigma = C.sigmalist[jM][1];
        dsigma = C.sigmalist[jM][2];
        
        outfile << M << "   " << sigma << "   " << dsigma << endl;
    }
    outfile.close();
    outfile.open("HMF_CDM.dat");
    for (int jz = 0; jz < C.Nz; jz++) {
        z = C.zlist[jz];
        for (int jM = 0; jM < C.NM; jM++) {
            M = C.HMFlist[jz][jM][1];
            HMF = C.HMFlist[jz][jM][2];
            dotM = C.dotMlist[jz][jM][2];
            
            outfile << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << endl;
        }
    }
    outfile.close();
    
    // output the FDM halo mass functions
    if (C.m22list.size() > 0) {
        outfile.open("sigma_FDM.dat");
        for (int jm = 0; jm < C.m22list.size(); jm++) {
            m22 = C.m22list[jm];
            for (int jM = 0; jM < C.NM; jM++) {
                M = C.FDMsigmalist[jm][jM][0];
                sigma = C.FDMsigmalist[jm][jM][1];
                dsigma = C.FDMsigmalist[jm][jM][2];
                
                outfile << m22 << "   " << M << "   " << sigma << "   " << dsigma << endl;
            }
        }
        outfile.close();
        
        outfile.open("HMF_FDM.dat");
        for (int jm = 0; jm < C.m22list.size(); jm++) {
            m22 = C.m22list[jm];
            for (int jz = 0; jz < C.Nz; jz++) {
                z = C.zlist[jz];
                for (int jM = 0; jM < C.NM; jM++) {
                    M = C.FDMHMFlist[jm][jz][jM][1];
                    HMF = C.FDMHMFlist[jm][jz][jM][2];
                    dotM = C.FDMdotMlist[jm][jz][jM][2];
                    
                    outfile << m22 << "   " << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << endl;
                }
            }
        }
        outfile.close();
    }
    
    // output the WDM halo mass functions
    if (C.m3list.size() > 0) {
        outfile.open("sigma_WDM.dat");
        for (int jm = 0; jm < C.m3list.size(); jm++) {
            m3 = C.m3list[jm];
            for (int jM = 0; jM < C.NM; jM++) {
                M = C.WDMsigmalist[jm][jM][0];
                sigma = C.WDMsigmalist[jm][jM][1];
                dsigma = C.WDMsigmalist[jm][jM][2];
                
                outfile << m3 << "   " << M << "   " << sigma << "   " << dsigma << endl;
            }
        }
        outfile.close();
        
        outfile.open("HMF_WDM.dat");
        for (int jm = 0; jm < C.m3list.size(); jm++) {
            m3 = C.m3list[jm];
            for (int jz = 0; jz < C.Nz; jz++) {
                z = C.zlist[jz];
                for (int jM = 0; jM < C.NM; jM++) {
                    M = C.WDMHMFlist[jm][jz][jM][1];
                    HMF = C.WDMHMFlist[jm][jz][jM][2];
                    dotM = C.WDMdotMlist[jm][jz][jM][2];
                    
                    outfile << m3 << "   " << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << endl;
                }
            }
        }
        outfile.close();
    }
        
    cout << "Generating/reading lensing amplifications..." << endl;
    
    // number of bins in distribution of P^1(kappa)
    int Nkappa = 40;
    // number of realizations
    int Nreal = 100000000;
    // number of lnmu bins
    int Nbins = 200;
    // lensing source radius in kpc
    double rS = 10.0;
    // list of redshift values at which P(lnmu) is computed
    C.Zlist = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5};
    C.Plnmuz = getPlnmu(C, rS, Nkappa, Nreal, Nbins);
    
    // read UV luminosity data files, {MUV, Phi, +sigmaPhi, -sigmaPhi}
    vector<string> datafiles {"UVLF_2102.07775.txt", "UVLF_2108.01090.txt", "UVLF_2403.03171.txt"};
    vector<vector<double> > data = readUVdata(datafiles);
    
    cout << "Computing UV luminosity fit (" << data.size() << " data points)..." << endl;
    
    // max number of steps in each chain
    int Nsteps = 2400;
    // burn-in
    int Nbi = 1200;
    // number of chains
    int Nchains = 200;
    
    // flat priors, {M_c, epsilon, alpha, beta, gamma, z_break, log_10(m)}
    vector<vector<double> > priors;
    if (C.m22list.size() > 0) {
        priors = {{3.0e11, 8.0e11}, {0.015, 0.045}, {0.6, 1.2}, {0.1, 0.6}, {0.0, 4.0}, {9.0, 11.0}, {log10(C.m22list.front()), log10(C.m22list.back())}};
    } else {
        priors = {{3.0e11, 8.0e11}, {0.015, 0.045}, {0.6, 1.2}, {0.1, 0.6}, {0.0, 4.0}, {9.0, 11.0}, {log10(C.m3list.front()), log10(C.m3list.back())}};
    }
    
    // random walk step sizes
    vector<double> steps(priors.size(),0.0);
    for (int j = 0; j < steps.size(); j++) {
        steps[j] = (priors[j][1]-priors[j][0])/20.0;
    }
    
    // random number generator
    rgen mt(time(NULL));
    
    // output the MCMC chains and find the best fit
    outfile.open("MCMCchains.dat");
    int jmax = 0;
    double logLmax = 0.0;
    vector<double> initial(priors.size(),0.0);
    vector<double> bf = initial;
    vector<vector<double> > chain;
    for (int j = 0; j < Nchains; j++) {
        cout << j << endl;
        
        // generate initial point inside the prior ranges
        for (int jp = 0; jp <  initial.size(); jp++) {
            initial[jp] = randomreal(priors[jp][0],priors[jp][1]);
        }
        
        // generate an MCMC chain
        chain = mcmc_sampling(C, data, initial, steps, priors, Nsteps, Nbi, mt);
        
        // find the best fit
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
    
    //vector<double> bf {6.0e11, 0.035, 0.95, 0.31, 1.1, 10.0, 0.29};
    
    // output the UV luminosity function for the best fit
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, 1.0e9, bf[0], bf[1], bf[2], bf[3], bf[4], bf[5], pow(10.0,bf[6]));
    //vector<vector<vector<double> > > PhiUVlist = PhiUV(C, 1.0e9, bf[0], bf[1], bf[2], bf[3], bf[4], bf[5]);
    outfile.open("UVluminosity.dat");
    double MUV, Phi0, Phi1, Phi2;
    for (int jZ = 0; jZ < C.Zlist.size(); jZ++) {
        z = C.Zlist[jZ];
        for (int jM = 0; jM < C.NM; jM++) {
            MUV = PhiUVlist[jZ][jM][0];
            M = PhiUVlist[jZ][jM][1];
            Phi0 = PhiUVlist[jZ][jM][2]; // no dust + no lensing
            Phi1 = PhiUVlist[jZ][jM][3]; // dust + no lensing
            Phi2 = PhiUVlist[jZ][jM][4]; // dust + lensing
            
            outfile << z << "   " << MUV << "   " << M << "   " << max(1.0e-64,Phi0) << "   " << max(1.0e-64,Phi1) << "   " << max(1.0e-64,Phi2) << endl;
        }
    }
    outfile.close();

    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
