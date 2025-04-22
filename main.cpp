#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(6) << fixed;
    
    const int doUVfit = atoi(argv[1]);
    const int dm = atoi(argv[2]); // 0: CDM, 1: FDM, 2: WDM
    
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
    C.Mmin = 1.0e6;
    C.Mmax = 1.0e17;
    C.NM = 2000;
    
    // redshifts
    C.zmin = 0.01;
    C.zmax = 36.01;
    C.Nz = 180;
    
    if (dm == 1) {
        // FDM masses in 10^-22 eV
        C.m22min = 1.0;
        C.m22max = 2000.0;
        C.Nm22 = 50;
    }
    if (dm == 2) {
        // WDM masses in keV
        C.m3min = 0.7;
        C.m3max = 40.0;
        C.Nm3 = 50;
    }
    if (dm == 3) {
        // kc values in 1/kpc
        C.kcmin = 0.006;
        C.kcmax = 2.0;
        C.Nkc = 50;
    }
    
    C.initialize();
    ofstream outfile;
   
    cout << "Computing halo mass functions..." << endl;
    
    double m22, m3, kc, z, M, sigma, dsigma, HMF, dotM;
    
    if (dm == 0) {
        C.CDM_halos();
        
        outfile.open("sigma_CDM.dat");
        for (int jM = 0; jM < C.NM; jM+=10) {
            M = C.sigmalist[jM][0];
            sigma = C.sigmalist[jM][1];
            dsigma = C.sigmalist[jM][2];
            
            outfile << M << "   " << sigma << "   " << dsigma << endl;
        }
        outfile.close();
        
        outfile.open("HMF_CDM.dat");
        for (int jz = 0; jz < C.Nz; jz++) {
            z = C.zlist[jz];
            for (int jM = 0; jM < C.NM; jM+=10) {
                M = C.HMFlist[jz][jM][1];
                HMF = C.HMFlist[jz][jM][2];
                dotM = C.HMFlist[jz][jM][3];
                
                outfile << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << endl;
            }
        }
        outfile.close();
    }
        
    if (dm == 1) {
        C.FDM_halos();
        
        outfile.open("sigma_FDM.dat");
        for (int jm = 0; jm < C.m22list.size(); jm++) {
            m22 = C.m22list[jm];
            for (int jM = 0; jM < C.NM; jM+=10) {
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
                for (int jM = 0; jM < C.NM; jM+=10) {
                    M = C.FDMHMFlist[jm][jz][jM][1];
                    HMF = C.FDMHMFlist[jm][jz][jM][2];
                    dotM = C.FDMHMFlist[jm][jz][jM][3];
                    
                    outfile << m22 << "   " << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << endl;
                }
            }
        }
        outfile.close();
    }
    
    if (dm == 2) {
        C.WDM_halos();
        
        outfile.open("sigma_WDM.dat");
        for (int jm = 0; jm < C.m3list.size(); jm++) {
            m3 = C.m3list[jm];
            for (int jM = 0; jM < C.NM; jM+=10) {
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
                for (int jM = 0; jM < C.NM; jM+=10) {
                    M = C.WDMHMFlist[jm][jz][jM][1];
                    HMF = C.WDMHMFlist[jm][jz][jM][2];
                    dotM = C.WDMHMFlist[jm][jz][jM][3];
                    
                    outfile << m3 << "   " << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << endl;
                }
            }
        }
        outfile.close();
    }
    
    if (dm == 3) {
        C.EDM_halos();
        
        outfile.open("sigma_EDM.dat");
        for (int jm = 0; jm < C.kclist.size(); jm++) {
            kc = C.kclist[jm];
            for (int jM = 0; jM < C.NM; jM+=10) {
                M = C.EDMsigmalist[jm][jM][0];
                sigma = C.EDMsigmalist[jm][jM][1];
                dsigma = C.EDMsigmalist[jm][jM][2];
                
                outfile << kc << "   " << M << "   " << sigma << "   " << dsigma << endl;
            }
        }
        outfile.close();
        
        outfile.open("HMF_EDM.dat");
        for (int jm = 0; jm < C.kclist.size(); jm++) {
            kc = C.kclist[jm];
            for (int jz = 0; jz < C.Nz; jz++) {
                z = C.zlist[jz];
                for (int jM = 0; jM < C.NM; jM+=10) {
                    M = C.EDMHMFlist[jm][jz][jM][1];
                    HMF = C.EDMHMFlist[jm][jz][jM][2];
                    dotM = C.EDMHMFlist[jm][jz][jM][3];
                    
                    outfile << kc << "   " << z << "   " << M << "   " << max(1.0e-99,HMF) << "   " << dotM << endl;
                }
            }
        }
        outfile.close();
    }
    
    cout << "Generating/reading lensing amplifications..." << endl;
    
    // number of bins in distribution of P^1(kappa)
    int Nkappa = 50;
    // number of realizations
    int Nreal = 2e8;
    // number of lnmu bins
    int Nbins = 400;
    // lensing source radius in kpc
    double rS = 10.0;
    // list of redshift values at which dP(lnmu)/dlnmu is computed
    //C.Zlist = {1.0, 2.0, 4.0, 8.0, 16.0};
    C.Zlist = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.5, 14.5, 17.0, 25.0};
    C.Plnmuz = getPlnmu(C, rS, Nkappa, Nreal, Nbins);
    
    // parameters: (logMt, Mc, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0, sigmaUV, logm)
    vector<double> bf = {7.89, 3.68e11, 0.0616, 0.915, 0.398, 0.222, 10.75, 0.353, 11.41, 25.25, 0.058, 0.7};
    //vector<double> bf = {6.02, 3.94e11, 0.0620, 0.737, 0.452, 0.201, 10.77, 0.355, 10.44, 24.74, 0.060, 0.08};
    //vector<double> bf = {7.89, 3.68e11, 0.0616, 0.915, 0.398, 0.222, 10.75, 0.353, 11.41, 25.25, 0.058, -1.6};
    
    string filename;
    if (doUVfit == 1) {
        cout << "Computing UV luminosity fit..." << endl;
        
        if (dm == 0) {
            filename = "MCMCchains_CDM.dat";
        }
        if (dm == 1) {
            filename = "MCMCchains_FDM.dat";
        }
        if (dm == 2) {
            filename = "MCMCchains_WDM.dat";
        }
        if (dm == 3) {
            filename = "MCMCchains_EDM.dat";
        }
        
        // max number of steps in each chain
        int Nsteps = 2000;
        // burn-in
        int Nbi = 1000;
        // number of chains
        int Nchains = 200;
        // flat priors
        vector<vector<double> > priors;
        if (C.m22list.size() > 0) {
            priors = {{6.0,10.0}, {3.0e11, 5.0e11}, {0.0578, 0.0658}, {0.6, 1.1}, {0.2, 0.6}, {0.05, 0.6}, {9.8, 11.4}, {0.05, 0.7}, {7.0, 25.0}, {20.0, 36.0}, {0.05, 0.2}, {log10(C.m22list.front()), log10(C.m22list.back())}};
        } else if (C.m3list.size() > 0) {
            priors = {{6.0,10.0}, {3.0e11, 5.0e11}, {0.0578, 0.0658}, {0.6, 1.1}, {0.2, 0.6}, {0.05, 0.6}, {9.8, 11.4}, {0.05, 0.7}, {7.0, 25.0}, {20.0, 36.0}, {0.05, 0.2}, {log10(C.m3list.front()), log10(C.m3list.back())}};
        } else if (C.kclist.size() > 0) {
            priors = {{6.0,10.0}, {3.0e11, 5.0e11}, {0.0578, 0.0658}, {0.6, 1.1}, {0.2, 0.6}, {0.05, 0.6}, {9.8, 11.4}, {0.05, 0.7}, {7.0, 25.0}, {20.0, 36.0}, {0.05, 0.2}, {log10(C.kclist.front()), log10(C.kclist.back())}};
        } else {
            priors = {{6.0,10.0}, {3.0e11, 5.0e11}, {0.0578, 0.0658}, {0.6, 1.1}, {0.2, 0.6}, {0.05, 0.6}, {9.8, 11.4}, {0.05, 0.7}, {7.0, 25.0}, {20.0, 36.0}, {0.05, 0.2}};
        }
        
        // random walk step sizes
        vector<double> steps(priors.size(),0.0);
        for (int j = 0; j < steps.size(); j++) {
            steps[j] = (priors[j][1]-priors[j][0])/14.0;
        }
        
        bf = UFfit(C, priors, steps, Nsteps, Nbi, Nchains, filename);
    }
    
    // output the UV luminosity function for the best fit
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, bf[0], bf[1], bf[2], bf[3], bf[4], bf[5], bf[6], bf[7], bf[8], bf[9], bf[10], bf[11]);
    cout << loglikelihood(C,bf) << endl;
    
    if (dm == 0) {
        filename = "UVluminosity_CDM.dat";
    }
    if (dm == 1) {
        filename = "UVluminosity_FDM.dat";
    }
    if (dm == 2) {
        filename = "UVluminosity_WDM.dat";
    }
    if (dm == 3) {
        filename = "UVluminosity_EDM.dat";
    }
    
    outfile.open(filename);
    outfile << scientific << setprecision(12);

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
