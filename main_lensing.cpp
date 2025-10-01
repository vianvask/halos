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
    rgen mt(time(NULL)); // random number generator

    const int dm = atoi(argv[1]); // 0: CDM, 1: FDM, 2: WDM, 3: EDM
    const int lensing = atoi(argv[2]); // 0: without lensing, 1: with lensing
    const int catl = atoi(argv[3]); // 0: DES-5yr SN catalogue, 1: SMBH catalogue, 2: NS catalogue, 3: combined SMBH+NS catalogue
    const double fstep = atof(argv[4]); // stepsize relative to the prior range

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
    C.zmax = 10.01;
    C.Nz = 40;
    
    // initialize benchmar CDM halo mass function
    C.initialize(0);
    
    ofstream outfile, outfile1, outfile2;
    
    int Nreal = 2e4; // realizations
    int Nhalos = 40; // number of halos in each realization
    int Nbins = 12; // P(lnmu) bins
    
    //Plnmuf(C, 1.0, Nhalos, Nreal, Nbins, mt, 1);
    
    double z, lnmu, DL0, DL, sigmaDL;
    vector<vector<double> > Plnmu;
    if (!fileExists("Plnmuz.dat")) {
        cout << "Generating lensing amplifications..." << endl;
        
        outfile1.open("Plnmuz.dat");
        outfile2.open("PDL.dat");
        
        C.Zlist = {0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

        vector<vector<double> > PDL;
        vector<double> sample, DLCL(2,0.0);
        for (int jz = 0; jz < C.Zlist.size(); jz++) {
            z = C.Zlist[jz];
            cout << "z = " << z << endl;
            
            // compute and output the distribution of lnmu
            Plnmu = Plnmuf(C, z, Nhalos, Nreal, Nbins, mt, 0);
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
    
    if (!fileExists("catalogueSMBH.dat")) {
        cout << "Generating SMBHB catalogue..." << endl;
        
        outfile.open("catalogueSMBH.dat");
        
        int Nheavy = 10; // size of the LISA SMBH catalogue with EM counterparts
        double zthr = 10.0; // threshold z
        double FsigmaDL = 0.003; // relative error in D_L for LISA SMBH binaries
        
        vector<double> zSMBHlist = readdata("zSMBHlist.dat"); // catalogue of LISA SMBH redshifts
        shuffle(zSMBHlist.begin(), zSMBHlist.end(), mt);
        zSMBHlist.erase(zSMBHlist.begin() + 10, zSMBHlist.end());
        
        int jD = 0, jz = 0;
        while (jD < Nheavy) {
            // generate redshift
            z = zSMBHlist[jz];
            jz++;
            
            if (z < zthr) {
                // generate amplification
                Plnmu = Plnmuf(C, z, Nhalos, Nreal, Nbins, mt, 0);
                lnmu = sampleFromPDF(Plnmu, 1, mt)[0];
                
                // lensed luminosity distance
                DL0 = C.DL(z)/exp(lnmu/2.0);
                
                // observed luminosity distance
                sigmaDL = FsigmaDL*DL0;
                DL = normal_distribution<double>(DL0, sigmaDL)(mt);
                                
                if (DL > 0.0) {
                    jD++;
                    cout << jD << endl;
                    outfile << z << "   " << DL << "   " << sigmaDL << endl;
                }
            }
        }
        outfile.close();
    }
    
    if (!fileExists("catalogueNS.dat")) {
        cout << "Generating NSB catalogue..." << endl;
        
        outfile.open("catalogueNS.dat");
        
        int Nlight = 1000; // size of the ET NS catalogue
        double zthr = 2.0; // threshold z
        double FsigmaDL = 0.03; // relative error in D_L for ET NS binaries
        double zeta = 2.6; // z dependence of NS binary mergers, N propto (1+z)^zeta
        
        // z distribution of the NS sources
        vector<vector<double> > CDFz(1000, vector<double> (2,0.0));
        CDFz[0][0] = C.zmin;
        CDFz[0][1] = 0.0;
        for (int jz = 1; jz < CDFz.size(); jz++) {
            z = pow(10.0, log10(C.zmin) + jz*(log10(zthr)-log10(C.zmin))/(1.0*CDFz.size()-1));
            CDFz[jz][0] = z;
            CDFz[jz][1] = CDFz[jz-1][1] + pz(C,z,zeta)*(CDFz[jz][0]-CDFz[jz-1][0]);
        }
        
        int jD = 0;
        while (jD < Nlight) {
            // generate redshift
            z = sampleFromCDF(CDFz, 1, mt)[0];
            
            // generate amplification
            Plnmu = Plnmuf(C, z, Nhalos, Nreal, Nbins, mt, 0);
            lnmu = sampleFromPDF(Plnmu, 1, mt)[0];
            
            // lensed luminosity distance
            DL0 = C.DL(z)/exp(lnmu/2.0);
            
            // observed luminosity distance
            sigmaDL = FsigmaDL*DL0;
            DL = normal_distribution<double>(DL0, sigmaDL)(mt);
                        
            if (DL > 0.0) {
                jD++;
                cout << jD << endl;
                outfile << z << "   " << DL << "   " << sigmaDL << endl;
            }
        }
        outfile.close();
    }
    
    cout << "Performing MCMC scan..." << endl;

    // MCMC scan parameters
    vector<double> par = {C.OmegaM, C.sigma8, C.h, 1.0};
    vector<vector<double> > priors = {{0.21, 0.42}, {0.5, 1.1}, {0.64, 0.71}, {-0.6, 1.4}};
    int Npar = par.size();
    
    // read catalogue
    vector<vector<double> > catalogue, catalogue2;
    vector<double> minmaxz, minmaxz2, Zlist2;
    if (catl < 3) {
        if (catl == 0) {
            catalogue = readdata("SNcatalogue.dat", 3);
            par = {0.33, 0.8, 0.7, 0.0};
            priors = {{0.1, 0.6}, {0.4, 1.4}, {0.699, 0.701}, {-0.1, 0.1}};
        } else if (catl == 1) {
            catalogue = readdata("catalogueSMBH.dat", 3);
        } else {
            catalogue = readdata("catalogueNS.dat", 3);
        }
        minmaxz = findMinMax(catalogue, 0);
        C.Zlist = linlist(minmaxz[0], minmaxz[1], 6);
    }
    else {
        catalogue = readdata("catalogueNS.dat", 3);
        
        minmaxz = findMinMax(catalogue, 0);
        C.Zlist = linlist(minmaxz[0], minmaxz[1], 6);
        
        catalogue2 = readdata("catalogueSMBH.dat", 3);
        catalogue.insert(catalogue.end(), catalogue2.begin(), catalogue2.end());
        
        minmaxz2 = findMinMax(catalogue2, 0);
        Zlist2 = linlist(minmaxz[1], minmaxz2[1], 4);
        C.Zlist.insert(C.Zlist.end(), Zlist2.begin(), Zlist2.end());
    }
        
    vector<double> steps(Npar, 0.0);
    for (int j = 0; j < Npar; j++) {
        steps[j] = fstep*(priors[j][1]-priors[j][0]);
    }
    
    // generate MCMC chains
    int Nsteps = 2000;
    int Nburnin = 200;
    for (int j = 0; j < 8; j++) {
        Hubble_diagram_fit(C, 3.0e8, catalogue, par, steps, priors, Nsteps, Nburnin, lensing, dm, mt, "lensing_chains_" + to_string(dm) + "_" + to_string(lensing) + "_" + to_string(catl) + "_c_" + to_string(j) + ".dat");
    }
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
