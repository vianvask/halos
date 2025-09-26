#include "cosmology.h"
#include "lensing.h"

// units: masses in solar masses, time in Myr, length in kpc

// z distribution of sources in the mock catalogue
double pz(cosmology &C, double z, double zeta) {
    return pow(1.0+z,zeta)*C.DVc(z);
}

int main (int argc, char *argv[]) {
    clock_t time_req = clock(); // timing
    cout << setprecision(8) << fixed;

    const int dm = atoi(argv[1]); // 0: CDM, 1: FDM, 2: WDM, 3: EDM
    const int lensing = atoi(argv[2]); // 1: without lensing, 2: with lensing
    const int threshold = atoi(argv[3]); // 1: threshold on D_L, 2: threshold on z, >2: DES-5yr SN catalogue, no threshold
    
    // z dependence of the source number, N propto (1+z)^zeta
    double zeta = 2.6;
    
    // benchmark cosmological parameters: PDG values
    cosmology C;
    C.OmegaM = 0.315;
    C.OmegaB = 0.0493;
    C.zeq = 3402.0;
    C.sigma8 = 0.811;
    C.h = 0.674;
    C.T0 = 2.7255;
    C.ns = 0.965;
    
    /*
    // comparison with T11
    C.OmegaM = 0.274;
    C.OmegaB = 0.046;
    C.zeq = 3402.0;
    C.sigma8 = 0.812;
    C.h = 0.705;
    C.T0 = 2.7255;
    C.ns = 0.96;
    */
    
    // halo masses
    C.Mmin = 1.0e7;
    C.Mmax = 1.0e17;
    C.NM = 100;
    
    // redshifts
    C.zmin = 0.01;
    C.zmax = 3.41;
    C.Nz = 34;
    
    // initialize benchmar CDM halo mass function
    C.initialize(0);
    
    ofstream outfile1, outfile2;
    rgen mt(time(NULL)); // random number generator
        
    int Nkappa = 40; // kappa bins in P^(1)(kappa)
    int Nreal = 2e4; // realizations
    int Nhalos = 10; // number of halos in each realization
    int Nbins = 12; // P(lnmu) bins
    double rS = 0.0; // lensing source radius in kpc
    
    // Plnmuf(C, 1.0, rS, Nhalos, Nreal, Nbins, mt, 1);
    // PlnmufA(C, 1.0, rS, Nkappa, Nreal, Nbins, mt, 1);
    
    double DLthr = 20.0e6; // detectability threshold, D_L^meas < 20 Gpc
    double zthr = 2.5; // detectability threshold, z_L^meas < 2.4
    
    double FsigmaDL = 0.03; // relative error in D_L
    int Nsources = 3000; // size of the catalogue
    
    double z, lnmu, DL0, DL, sigmaDL;
    vector<vector<double> > Plnmu;
    if (!fileExists("Plnmuz.dat")) {
        cout << "Generating lensing amplifications..." << endl;
        
        outfile1.open("Plnmuz.dat");
        outfile2.open("PDL.dat");
        
        C.Zlist = {0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4};

        vector<vector<double> > PDL;
        vector<double> sample, DLCL(2,0.0);
        for (int jz = 0; jz < C.Zlist.size(); jz++) {
            z = C.Zlist[jz];
            cout << "z = " << z << endl;
            
            // compute and output the distribution of lnmu
            Plnmu = Plnmuf(C, z, rS, Nhalos, Nreal, Nbins, mt, 0);
            for (int jb = 0; jb < Plnmu.size(); jb++) {
                outfile1 << z << "   " << Plnmu[jb][0] << "   " << Plnmu[jb][1] << endl;
            }
            
            // generate a sample of apparent D_L
            sample = sampleFromPDF(Plnmu, Nreal, mt);
            DL0 = C.DL(z);
            for (int j = 0; j < Nreal; j++) {
                sample[j] = DL0/exp(sample[j]/2.0);
            }
            
            // compute and output the variance and the 99% confidence interval of apparent D_L
            DLCL = confidenceInterval(sample, 0.99);
            outfile2 << z << "   " << DL0 << "   " << variance(sample) << "   " << DLCL[1] << "   " << DLCL[2] << endl;
        }
        outfile1.close();
        outfile2.close();
    }
    
    if (!fileExists("catalogueD.dat") || !fileExists("catalogueZ.dat")) {
        cout << "Generating source catalogues..." << endl;
        
        // z distribution of the sources
        vector<vector<double> > CDFz(1000, vector<double> (2,0.0));
        CDFz[0][0] = C.zmin;
        CDFz[0][1] = 0.0;
        for (int jz = 1; jz < CDFz.size(); jz++) {
            z = pow(10.0, log10(C.zmin) + jz*(log10(C.zmax)-log10(C.zmin))/(1.0*CDFz.size()-1));
            CDFz[jz][0] = z;
            CDFz[jz][1] = CDFz[jz-1][1] + pz(C,z,zeta)*(CDFz[jz][0]-CDFz[jz-1][0]);
        }
        
        outfile1.open("catalogueD.dat");
        outfile2.open("catalogueZ.dat");
        
        // generate a catalogue of sources
        int jD = 0, jZ = 0;
        while (jD < Nsources || jZ < Nsources) {
            // generate redshift
            z = sampleFromCDF(CDFz, 1, mt)[0];
            
            // generate amplification
            Plnmu = Plnmuf(C, z, rS, Nhalos, Nreal, Nbins, mt, 0);
            lnmu = sampleFromPDF(Plnmu, 1, mt)[0];
            
            // lensed luminosity distance
            DL0 = C.DL(z)/exp(lnmu/2.0);
            
            // observed luminosity distance
            sigmaDL = FsigmaDL*DL0;
            DL = normal_distribution<double>(DL0, sigmaDL)(mt);
            
            if (DL > 0.0) {
                if (DL < DLthr && jD < Nsources) {
                    jD++;
                    cout << jD << "   " << jZ << endl;
                    outfile1 << z << "   " << DL << "   " << sigmaDL << endl;
                }
                
                if (z < zthr && jZ < Nsources) {
                    jZ++;
                    cout << jD << "   " << jZ << endl;
                    outfile2 << z << "   " << DL << "   " << sigmaDL << endl;
                }
            }
        }
        outfile1.close();
        outfile2.close();
    }
    
    vector<vector<double> > catalogue;
    if (threshold == 1) {
        catalogue = readdata("catalogueD.dat", 3);
    } else if (threshold == 2) {
        catalogue = readdata("catalogueZ.dat", 3);
        DLthr = 2.0*DLthr;
    } else {
        catalogue = readdata("SNcatalogue.dat", 3);
        vector<double> minmaxD = findMinMax(catalogue, 1);
        DLthr = 2.0*minmaxD[1];
        cout << DLthr << endl;
    }
    
    // make a new list of z values
    vector<double> minmaxz = findMinMax(catalogue, 0);
    cout << minmaxz[0] << "   " << minmaxz[1] << endl;
    C.Zlist = linlist(minmaxz[0], minmaxz[1], 6);
        
    cout << "Performing MCMC scan..." << endl;
    
    // MCMC scan over Omega_M, sigma_8 and h:
    vector<double> par = {C.OmegaM, C.sigma8, C.h, 1.0};
    //vector<double> par = {0.33, C.sigma8, 0.7, 0.0};
    int Npar = par.size();
    
    vector<vector<double> > priors = {{0.28, 0.38}, {0.65, 1.05}, {0.646, 0.69}, {-0.6, 1.4}};
    //vector<vector<double> > priors = {{0.1, 0.6}, {0.4, 1.4}, {0.699, 0.701}, {-0.1, 0.1}};
    vector<double> steps(Npar,0.0);
    for (int j = 0; j < Npar; j++) {
        steps[j] = (priors[j][1]-priors[j][0])/16.0;
    }
    int Nsteps = 2000;
    int Nburnin = 200;
    
    outfile1.open("lensing_chains_" + to_string(dm) + "_" + to_string(lensing) + "_" + to_string(threshold) + ".dat");
    vector<vector<double> > chain;
    for (int jc = 0; jc < 8; jc++) {
        chain = Hubble_diagram_fit(C, DLthr, catalogue, par, steps, priors, Nsteps, Nburnin, lensing, dm, mt);
        
        // output the chain
        for (int js = 0; js < Nsteps; js++) {
            for (int jp = 0; jp < Npar+1; jp++) {
                outfile1 << chain[js][jp] << "   ";
            }
            outfile1 << endl;
        }
    }
    outfile1.close();
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
