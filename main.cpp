#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(6) << fixed;
    
    const int doUVfit = atoi(argv[1]);
    const int NmF = atoi(argv[2]);
    const int NmW = atoi(argv[3]);
    
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
    C.zmax = 40.01;
    C.Nz = 160;
    
    // FDM masses in 10^-22 eV
    C.m22min = 0.8;
    C.m22max = 80.0;
    C.Nm22 = NmF;
    
    // WDM masses in keV
    C.m3min = 0.5;
    C.m3max = 50.0;
    C.Nm3 = NmW;
    
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
    vector<double> bf = {8.0, 4.1e11, 0.063, 0.94, 0.39, 0.17, 10.8, 0.33, 10.3, 30.4, 0.056, 0.0};
    
    if (doUVfit == 1) {
        cout << "Computing UV luminosity fit..." << endl;
        
        // max number of steps in each chain
        int Nsteps = 2000;
        // burn-in
        int Nbi = 1000;
        // number of chains
        int Nchains = 100;
        // flat priors
        vector<vector<double> > priors;
        if (C.m22list.size() > 0) {
            priors = {{6.0,9.0}, {3.2e11, 4.6e11}, {0.058, 0.068}, {0.76 , 1.1}, {0.24, 0.6}, {0.05, 0.7}, {9.5, 11.5}, {0.05,  1.0}, {8.0, 25.0}, {15.0, 40.0}, {0.05, 0.25}, {log10(C.m22list.front()), log10(C.m22list.back())}};
        } else if (C.m3list.size() > 0) {
            priors = {{6.0,9.0}, {3.2e11, 4.6e11}, {0.058, 0.068}, {0.76 , 1.1}, {0.24, 0.6}, {0.05, 0.7}, {9.5, 11.5}, {0.05,  1.0}, {8.0, 25.0}, {15.0, 40.0}, {0.05, 0.25}, {log10(C.m3list.front()), log10(C.m3list.back())}};
        } else {
            priors = {{6.0,9.0}, {3.2e11, 4.6e11}, {0.058, 0.068}, {0.76 , 1.1}, {0.24, 0.6}, {0.05, 0.7}, {9.5, 11.5}, {0.05,  1.0}, {8.0, 25.0}, {15.0, 40.0}, {0.05, 0.25}};
        }
        
        // random walk step sizes
        vector<double> steps(priors.size(),0.0);
        for (int j = 0; j < steps.size(); j++) {
            steps[j] = (priors[j][1]-priors[j][0])/14.0;
        }
        
        bf = UFfit(C, priors, steps, Nsteps, Nbi, Nchains);
    }
    
    // output the UV luminosity function for the best fit
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, bf[0], bf[1], bf[2], bf[3], bf[4], bf[5], bf[6], bf[7], bf[8], bf[9], bf[10], bf[11]);
    cout << loglikelihood(C,bf) << endl;
    
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
    
    /*
    cout << "Evolving stellar masses and BH masses..." << endl;
    
    outfile.open("BHmassStellarmass.dat");
    vector<vector<double> > MBHlist(C.Nz, vector<double> (C.NM,0.0));
    vector<vector<double> > Mstlist = C.evolvestellarmass(bf[0], bf[1], bf[2], bf[3], bf[4]);
    
    double MBH, Mst;
    for (int jz = 0; jz < C.Nz; jz++) {
        z = C.zlist[jz];
        for (int jM = 0; jM < C.NM; jM++) {
            M = C.sigmalist[jM][0];
            MBH = MBHlist[jz][jM];
            Mst = Mstlist[jz][jM];
                        
            outfile << z << "   " << M << "   " << MBH << "   " << Mst << endl;
        }
    }
    outfile.close();
    */
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
