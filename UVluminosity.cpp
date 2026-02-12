#include "cosmology.h"
#include "UVluminosity.h"

/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                       lensing magnifications                                                                   */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// read lensing amplification distribution
vector<vector<vector<double> > > getPlnmu(fs::path filename) {
    vector<vector<double> > Plnmu;
    vector<vector<vector<double> > > Plnmuz;
    
    double z;
    vector<double> tmp(2, 0.0);
    
    ifstream infile;
    infile.open(filename);
    if (infile) {
        int jA = 0;
        double A;
        z = 0;
        while (infile >> A) {
            if (jA == 0) {
                if (A > z) {
                    if (z > 0) {
                        Plnmuz.push_back(Plnmu);
                        Plnmu.clear();
                    }
                    z = A;
                }
                jA++;
            } else {
                tmp[jA-1] = A;
                jA++;
            }
            if (jA == 3) {
                Plnmu.push_back(tmp);
                jA = 0;
            }
        }
        Plnmuz.push_back(Plnmu);
        Plnmu.clear();
    }
    else {
        cout << "Plnmu.dat missing." << endl;
    }
    infile.close();
        
    double dlnmu;
    for (int jz = 0; jz < Plnmuz.size(); jz++) {
        
        // find the element where the second component is 0
        auto vec = Plnmuz[jz];
        auto it = std::find_if(vec.begin(), vec.end(), [](const std::vector<double> &v) { return v.size() > 1 && v[1] == 0.0; });
        
        // erase from that point to the end
        vec.erase(it, vec.end());
        
        // extrapolate
        dlnmu = vec[vec.size()-1][0] - vec[vec.size()-2][0];
        while (vec.back()[0] < 1000.0) {
            tmp[0] = vec.back()[0] + dlnmu;
            tmp[1] = vec.back()[1]*exp(-3.0*dlnmu);
            vec.push_back(tmp);
        }
        Plnmuz[jz] = vec;
    }
    return Plnmuz;
}



/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                              UVLF model                                                                        */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// UV luminosity conversion factor
double kappaUV(double z, double gamma, double zc, double fkappa) {
    double kappaUV0 = 1.15e-22; // Msun s erg^-1 Myr^-1
    return kappaUV0*((1.0+fkappa)-(1.0-fkappa)*tanh((z-zc)/gamma))/2.0;
}

// change alpha and beta at z>ze
double enh(double z, double ze, double z0) {
    if (z0 > z && z > ze) {
        return (z0-z)/(z0-ze);
    }
    if (z >= z0) {
        return 0.0;
    }
    return 1.0;
}

// the UV magnitude
double MUV(cosmology &C, double z, double M, double dotM, double Mc, double Mt, double epsilon, double alpha, double beta, double gamma, double zc, double fkappa, double ze, double z0) {
    double alphaz = alpha*enh(z, ze, z0);
    double betaz = beta*enh(z, ze, z0);
    return 51.63 - 1.08574*log(C.fstar(z,M,Mc,Mt,epsilon,alphaz,betaz)*max(1.0e-99,dotM)/kappaUV(z,gamma,zc,fkappa));
}

// derivative of the UV magnitude, dM_UV/dM
double DMUV(cosmology &C, double z, double M, double dotM, double DdotM, double Mc, double Mt, double epsilon, double alpha, double beta) {
    return -1.08574*(DdotM/max(1.0e-99,dotM) + C.Dfstarperfstar(z,M,Mc,Mt,epsilon,alpha,beta));
}

// dust extinction, MUV -> MUV - AUV (see 1406.1503 and Table 3 of 1306.2950)
double AUV(double MUV, double z) {
    double C0 = 4.4, C1 = 2.0, sigmabeta = 0.34;
    double beta = -exp(0.17*(19.5+MUV)/(1.54+0.075*z))*(1.54+0.075*z);
    double AUV0 = C0 + 0.2*log(10)*pow(C1*sigmabeta,2.0) + C1*beta; // luminosity based average (1904.07238)
    double s = 3.0; // approach smoothly AUV = 0 (compare with Fig 1 of  1904.07238)
    return log(exp(s*AUV0) + 1.0)/s;
}

// UV luminosity function as a function of halo mass M, dP(MUV|M)/dMUV = delta(MUV - MUV(M)) (see e.g. 1906.06296)
double UVLF(cosmology &C, double z, double M, double dotM, double DdotM, double dndlnM, double Mc, double Mt, double epsilon, double alpha, double beta, double ze, double z0) {
    double alphaz = alpha*enh(z, ze, z0);
    double betaz = beta*enh(z, ze, z0);
    return -dndlnM/M/DMUV(C,z,M,dotM,DdotM,Mc,Mt,epsilon,alphaz,betaz);
}

// distribution of emitted UV magnitudes
double pMUV(double MUV, double MUVbar, double sigmaUV) {
    return 1.0/(sqrt(2.0*PI)*sigmaUV)*exp(-pow(MUV-MUVbar,2.0)/(2.0*pow(sigmaUV,2.0)));
}

// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
vector<vector<vector<double> > > PhiUV(cosmology &C, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma, double zc, double fkappa, double ze, double z0, double sigmaUV) {
    
    vector<vector<vector<double> > > PhiUVlensed(C.Zlist.size(), vector<vector<double> > (C.NM, vector<double> (5,0.0)));
    vector<vector<double> > PhiUVlist(C.NM, vector<double> (4,0.0));
    
    double z, dlnmu, Mj, dotMj, DdotMj, dndlnMj, MUVj, UVLFj, AUVj, PhiUVj;
    double Mi, dotMi, dndlnMi, MUVi, Mi2, dotMi2, dndlnMi2, MUVi2;
    vector<vector<double> > Plnmu(C.Plnmuz[0].size(), vector<double> (2,0.0));
    
    int jzp;
    for (int jZ = 0; jZ < C.Zlist.size(); jZ++) {
        z = C.Zlist[jZ];
        Plnmu = C.Plnmuz[jZ];
        
        // find index of z in the longer list of z values
        jzp = lower_bound(C.zlist.begin(), C.zlist.end(), z)- C.zlist.begin();
        if (jzp > 0 && C.zlist[jzp]-z > z-C.zlist[jzp-1]) {
            jzp--;
        }
        
        // compute the unlensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            Mj = C.Mlist[jM];
            dndlnMj = C.HMFlist[jzp][jM][0];
            
            dotMj = C.HMFlist[jzp][jM][1];
            DdotMj = C.HMFlist[jzp][jM][2];
            
            MUVj = MUV(C, z, Mj, dotMj, Mc, Mt, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0);
            
            // TODO: check, there seems to be a problem with implementation of sigma_UV
            // integral over halo masses to account for the distribution of emitted UV magnitudes
            UVLFj = 0.0;
            if (sigmaUV > 0.0) {
                for (int iM = jM; iM < C.NM-1; iM++) { // values larger than the mean
                    Mi = C.Mlist[iM];
                    dndlnMi = C.HMFlist[jzp][iM][0];
                    dotMi = C.HMFlist[jzp][iM][1];
                    
                    Mi2 = C.Mlist[iM+1];
                    dndlnMi2 = C.HMFlist[jzp][iM+1][0];
                    dotMi2 = C.HMFlist[jzp][iM+1][1];
                    
                    MUVi = MUV(C, z, Mi, dotMi, Mc, Mt, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0);
                    MUVi2 = MUV(C, z, Mi2, dotMi2, Mc, Mt, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0);
                    
                    UVLFj += (dndlnMi*pMUV(MUVj,MUVi,sigmaUV) + dndlnMi2*pMUV(MUVj,MUVi2,sigmaUV))/2.0*(log(Mi2)-log(Mi));
                    
                    if (abs(MUVi-MUVj) > 5.0*sigmaUV) { // stop at 5sigma
                        iM = C.NM;
                    }
                }
                for (int iM = jM; iM > 0; iM--) { // values smaller than the mean
                    Mi = C.Mlist[iM];
                    dndlnMi = C.HMFlist[jzp][iM][0];
                    dotMi = C.HMFlist[jzp][iM][1];
                    
                    Mi2 = C.Mlist[iM-1];
                    dndlnMi2 = C.HMFlist[jzp][iM-1][0];
                    dotMi2 = C.HMFlist[jzp][iM-1][1];
                    
                    MUVi = MUV(C, z, Mi, dotMi, Mc, Mt, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0);
                    MUVi2 = MUV(C, z, Mi2, dotMi2, Mc, Mt, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0);
                    
                    UVLFj += (dndlnMi*pMUV(MUVj,MUVi,sigmaUV) + dndlnMi2*pMUV(MUVj,MUVi2,sigmaUV))/2.0*(log(Mi)-log(Mi2));
                    
                    if (abs(MUVi-MUVj) > 5.0*sigmaUV) { // stop at 5sigma
                        iM = 0;
                    }
                }
            }
            else { // dP(MUV|M)/dMUV = delta(MUV - MUV(M))
                UVLFj = UVLF(C, z, Mj, dotMj, DdotMj, dndlnMj, Mc, Mt, epsilon, alpha, beta, ze, z0);
            }
            
            PhiUVlist[jM][0] = MUVj;
            PhiUVlist[jM][1] = Mj;
            PhiUVlist[jM][2] = AUV(MUVj, z);
            PhiUVlist[jM][3] = UVLFj;
        }
        
        // compute the lensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            MUVj = PhiUVlist[jM][0];
            Mj = PhiUVlist[jM][1];
            AUVj = PhiUVlist[jM][2];
            
            // convolution
            PhiUVj = 0.0;
            dlnmu = Plnmu[1][0] - Plnmu[0][0];
            for (int jb = 0; jb < Plnmu.size(); jb++) {
                PhiUVj += interpolaten(MUVj - AUVj + 1.086*Plnmu[jb][0], PhiUVlist)[3]*Plnmu[jb][1]*dlnmu;
            }
            
            PhiUVlensed[jZ][jM][0] = MUVj;
            PhiUVlensed[jZ][jM][1] = Mj;
            PhiUVlensed[jZ][jM][2] = PhiUVlist[jM][3]; // no dust + no lensing
            PhiUVlensed[jZ][jM][3] = interpolaten(MUVj - AUVj, PhiUVlist)[3]; // dust + no lensing
            PhiUVlensed[jZ][jM][4] = PhiUVj; // dust + lensing
        }
    }
    return PhiUVlensed;
}


// generalization of the PhiUV function to different DM models
vector<vector<vector<double> > > PhiUV(cosmology &C, double logMt, double Mc, double epsilon, double alpha, double beta, double gamma, double zc, double fkappa, double ze, double z0, double sigmaUV, double logm, int dm) {
    
    double Mt = pow(10.0, logMt);
    double m = pow(10.0, logm);
    
    if (dm == 1) {
        int jm = lower_bound(C.m22list.begin(), C.m22list.end(), m)- C.m22list.begin();
        if (jm > 0 && C.m22list[jm]-m > m-C.m22list[jm-1]) {
            jm--;
        }
        C.HMFlist = C.FDMHMFlist[jm];
    }
    if (dm == 2) {
        int jm = lower_bound(C.m3list.begin(), C.m3list.end(), m)- C.m3list.begin();
        if (jm > 0 && C.m3list[jm]-m > m-C.m3list[jm-1]) {
            jm--;
        }
        C.HMFlist = C.WDMHMFlist[jm];
    }
    if (dm == 3) {
        int jm = lower_bound(C.kclist.begin(), C.kclist.end(), m)- C.kclist.begin();
        if (jm > 0 && C.kclist[jm]-m > m-C.kclist[jm-1]) {
            jm--;
        }
        C.HMFlist = C.EDMHMFlist[jm];
    }
    if (dm == 4) {
        int jm = lower_bound(C.Blist.begin(), C.Blist.end(), m)- C.Blist.begin();
        if (jm > 0 && C.Blist[jm]-m > m-C.Blist[jm-1]) {
            jm--;
        }
        C.HMFlist = C.BDMHMFlist[jm];
    }
    
    return PhiUV(C, Mt, Mc, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0, sigmaUV);
}


// write UV luminosity functions to a file for a benchmark case
void writeUVLF(cosmology &C, double logMt, double Mc, double epsilon, double alpha, double beta, double gamma, double zc, double fkappa, double ze, double z0, double sigmaUV, double logm, int dm) {
    
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, logMt, Mc, epsilon, alpha, beta, gamma, zc, fkappa, ze, z0, sigmaUV, logm, dm);
    
    string filename;
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
    if (dm == 4) {
        filename = "UVluminosity_B.dat";
    }
    
    ofstream outfile;
    outfile.open(C.outdir/filename);
    outfile << scientific << setprecision(12);

    double z, MUV, M, Phi0, Phi1, Phi2;
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
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                               UVLF fit                                                                         */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// read UV luminosity data, {z, MUV, Phi, +sigmaPhi, -sigmaPhi}
vector<vector<double> > readUVdata(vector<fs::path> filedir) {
    vector<vector<double> > data;
    vector<double> row(5,0.0);
    int jrow = 0;
    double A;
    
    ifstream infile;
    for (int j = 0; j < filedir.size(); j++) {
        infile.open(filedir[j]);
        if (infile) {
            while (infile >> A) {
                row[jrow] = A;
                jrow++;
                if (jrow == 5) {
                    jrow = 0;
                    data.push_back(row);
                }
            }
        } else {
            cout << "couldn't open " << filedir[j] << endl;
        }
        infile.close();
    }
    return data;
}

// linear uncertainties
double logNPDF2(double x, double mu, double sigmap, double sigmam) {
    if (x < mu) {
        return logNPDF(x,mu,sigmam) + log(2*sigmam/(sigmam+sigmap));
    }
    return logNPDF(x,mu,sigmap) + log(2*sigmap/(sigmam+sigmap));
}

// loglikelihood function
double loglikelihood(cosmology &C, vector<vector<double> > &data, vector<double> &params, int dm) {
        
    // compute the UV luminosity function
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10], params[11], dm);
    
    vector<vector<double> > PhiUVz;
    double z, MUV, meanPhi, sigmaPhip, sigmaPhim, Phi, logL = 0.0;
    for (int j = 0; j < data.size(); j++) {
        z = data[j][0];
        MUV = data[j][1];
        meanPhi = data[j][2];
        sigmaPhip = data[j][3];
        sigmaPhim = data[j][4];
        
        PhiUVz = PhiUVlist[lower_bound(C.Zlist.begin(), C.Zlist.end(), z) - C.Zlist.begin()];
        
        // Phi(M_UV), convert to Mpc^-3 as the data is in those units
        Phi = 1.0e9*interpolaten(MUV, PhiUVz)[4];
        
        logL += logNPDF2(Phi, meanPhi, sigmaPhip, sigmaPhim);
    }
    return logL;
}

// function to generate MCMC chains to fit the UV data
vector<double> UVLFfit(cosmology &C, vector<vector<double> > &priors, int Nsteps, int Nburnin, int Nchains, double xstep, int dm) {
    
    // random number generator
    rgen mt(time(NULL));
    
    // read UV luminosity data files, {MUV, Phi, +sigmaPhi, -sigmaPhi}
    vector<string> datafiles {"UVLF_2102.07775.txt", "UVLF_2108.01090.txt", "UVLF_2403.03171.txt", "UVLF_2503.15594.txt", "UVLF_2504.05893.txt", "UVLF_2505.11263.txt"};
    vector<fs::path> datadir(datafiles.size());
    for (int j = 0; j < datadir.size(); j++) {
        datadir[j] = C.outdir/datafiles[j];
    }
    vector<vector<double> > data = readUVdata(datadir);
    
    int Npar = priors.size();
    
    // random walk step sizes
    vector<double> steps(Npar,0.0);
    for (int j = 0; j < Npar; j++) {
        steps[j] = (priors[j][1] - priors[j][0])/xstep;
    }
    
    // output the MCMC chains and find the best fit
    vector<double> initial(Npar,0.0);
    vector<double> bf(Npar+1,0.0);
    vector<vector<double> > chain;
    for (int j = 0; j < Nchains; j++) {
        cout << j << endl;
        string filename = "UVLF_chains_" + to_string(dm) + "_c_" + to_string(j) + ".dat";
        
        // generate initial point inside the prior ranges
        for (int jp = 0; jp <  initial.size(); jp++) {
            initial[jp] = randomreal(priors[jp][0], priors[jp][1]);
        }
        
        // loglikelihood
        function<double(vector<double>&)> logpdf = [&C, &data, dm](vector<double> &par) {
            return loglikelihood(C, data, par, dm);
        };
        
        // no cut
        function<double(vector<double>&)> cut = [](vector<double> &par) {
            return 1.0;
        };
        
        // generate and output MCMC chain
        chain = MCMC_sampling(Nsteps, Nburnin, logpdf, initial, steps, priors, cut, mt, 1, 1, C.outdir/filename);
        
        // find the best fit
        for (int j = 0; j < chain.size(); j++) {
            if (chain[j].back() > bf.back()) {
                bf = chain[j];
            }
        }
    }
    
    return bf;
}
