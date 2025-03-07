#include "cosmology.h"
#include "UVluminosity.h"

double kappaUV = 1.15e-22; // Msun s erg^-1 Myr^-1

// star formation rate
double fstar(double M, double Mc, double epsilon, double alpha, double beta, double gamma, double z, double zbreak) {
    //return epsilon*max(1.0,pow(z/zbreak,gamma))*(alpha+beta)/(beta*pow(M/Mc,-alpha) + alpha*pow(M/Mc,beta));
    return epsilon*max(1.0,1.0+gamma*(z-zbreak))*(alpha+beta)/(beta*pow(M/Mc,-alpha) + alpha*pow(M/Mc,beta));
}

// derivative of the star formation rate, df_*/dM
double Dfstarperfstar(double M, double Mc, double epsilon, double alpha, double beta) {
    return beta*((alpha+beta)/(alpha*pow(M/Mc,alpha+beta)+beta) - 1)/M;
}

// the UV magnitude
double MUV(double M, double dotM, double Mc, double epsilon, double alpha, double beta, double gamma, double z, double zbreak) {
    return 51.63 - 1.08574*log(fstar(M,Mc,epsilon,alpha,beta,gamma,z,zbreak)*dotM/kappaUV);
}

// derivative of the UV magnitude, dM_UV/dM
double DMUV(double M, double dotM, double DdotM, double Mc, double epsilon, double alpha, double beta) {
    return -1.08574*(DdotM/dotM + Dfstarperfstar(M,Mc,epsilon,alpha,beta));
}

// dust extinction, MUV -> MUV - AUV (see 1406.1503 and Table 3 of 1306.2950)
double AUV(double MUV, double z) {
    double C0 = 4.4, C1 = 2.0, sigmabeta = 0.34;
    double beta = -exp(0.17*(19.5+MUV)/(1.54+0.075*z))*(1.54+0.075*z);
    return max(0.0, C0 + 0.2*log(10)*pow(C1*sigmabeta,2.0) + C1*beta);
}

// UV luminosity function as a function of halo mass M (see e.g. 1906.06296)
double UVLF(double M, double dotM, double DdotM, double dndlnM, double Mt, double Mc, double epsilon, double alpha, double beta) {
    return -exp(-Mt/M)*dndlnM/M/DMUV(M,dotM,DdotM,Mc,epsilon,alpha,beta);
}


// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
vector<vector<vector<double> > > PhiUV(cosmology &C, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma, double zbreak) {
    
    vector<vector<vector<double> > > PhiUVlensed(C.Zlist.size(), vector<vector<double> > (C.NM, vector<double> (5,0.0)));
    vector<vector<double> > PhiUVlist(C.NM, vector<double> (4,0.0));
    
    double z, dlnmu, Mj, dotMj, DdotMj, dndlnMj, MUVj, AUVj, PhiUVj;
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
            Mj = C.HMFlist[jzp][jM][1];
            dndlnMj = C.HMFlist[jzp][jM][2];
            
            dotMj = C.dotMlist[jzp][jM][2];
            DdotMj = C.dotMlist[jzp][jM][3];
            
            MUVj = MUV(Mj, dotMj, Mc, epsilon, alpha, beta, gamma, z, zbreak);
            
            PhiUVlist[jM][0] = MUVj;
            PhiUVlist[jM][1] = Mj;
            PhiUVlist[jM][2] = AUV(MUVj, z);
            PhiUVlist[jM][3] = UVLF(Mj, dotMj, DdotMj, dndlnMj, Mt, Mc, epsilon, alpha, beta);
        }
        
        // compute the lensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            MUVj = PhiUVlist[jM][0];
            Mj = PhiUVlist[jM][1];
            AUVj = PhiUVlist[jM][2];
            
            // convolution integral
            PhiUVj = 0.0;
            dlnmu = Plnmu[1][0]-Plnmu[0][0];
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

vector<vector<vector<double> > > PhiUV(cosmology &C, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma, double zbreak, double m) {
    
    // find index of z in the longer list of z values
    if (C.m22list.size() > 0) {
        int jm = lower_bound(C.m22list.begin(), C.m22list.end(), m)- C.m22list.begin();
        if (jm > 0 && C.m22list[jm]-m > m-C.m22list[jm-1]) {
            jm--;
        }
        C.HMFlist = C.FDMHMFlist[jm];
        C.dotMlist = C.FDMdotMlist[jm];
    } else {
        int jm = lower_bound(C.m3list.begin(), C.m3list.end(), m)- C.m3list.begin();
        if (jm > 0 && C.m3list[jm]-m > m-C.m3list[jm-1]) {
            jm--;
        }
        C.HMFlist = C.WDMHMFlist[jm];
        C.dotMlist = C.WDMdotMlist[jm];
    }
    
    return PhiUV(C, Mt, Mc, epsilon, alpha, beta, gamma, zbreak);
}

// read UV luminosity data, {z, MUV, Phi, +sigmaPhi, -sigmaPhi}
vector<vector<double> > readUVdata(vector<string> filenames) {
    vector<vector<double> > data;
    vector<double> row(5,0.0);
    int jrow = 0;
    double A;
    
    ifstream infile;
    for (int j = 0; j < filenames.size(); j++) {
        infile.open(filenames[j]);
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
            cout << "couldn't open " << filenames[j] << endl;
        }
        infile.close();
    }
    
    return data;
}


// normal distribution
double NPDF(double x, double mu, double sigma) {
    return exp(-pow((x-mu)/sigma,2.0)/2.0)/(sigma*sqrt(2.0*PI));
}
double NPDF2(double x, double mu, double sigmap, double sigmam) {
    if (x < mu) {
        return 2*sigmam*NPDF(x, mu, sigmam)/(sigmam+sigmap);
    }
    return 2*sigmap*NPDF(x, mu, sigmap)/(sigmam+sigmap);
}


// likelihood function
double loglikelihood(cosmology &C, vector<vector<double> > &data, vector<double> &params) {
        
    // compute the UV luminosity function
    vector<vector<vector<double> > > PhiUVlist;
    if (params.size() < 7) {
        PhiUVlist = PhiUV(C, 1.0e9, params[0], params[1], params[2], params[3], params[4], params[5]);
    } else {
        PhiUVlist = PhiUV(C, 1.0e9, params[0], params[1], params[2], params[3], params[4], params[5], pow(10.0,params[6]));
    }
    
    vector<vector<double> > PhiUVz;
    double z, MUV, meanPhi, sigmaPhip, sigmaPhim, Phi, logL = 0.0;
    for (int j = 0; j < data.size(); j++) {
        z = data[j][0];
        MUV = data[j][1];
        meanPhi = data[j][2];
        sigmaPhip = data[j][3];
        sigmaPhim = data[j][4];
        
        // z in Zlist
        PhiUVz = PhiUVlist[lower_bound(C.Zlist.begin(), C.Zlist.end(), z) - C.Zlist.begin()];
        
        // Phi(M_UV), convert to Mpc^-3 as data is in such units
        Phi = 1.0e9*interpolaten(MUV, PhiUVz)[4];
        
        logL += log(NPDF2(Phi, meanPhi, sigmaPhip, sigmaPhim));
    }
    return logL;
}


// uniform distributions as priors
double prior(double x, const vector<double> &bounds) {
    double lower = bounds[0], upper = bounds[1];
    if (x < lower || x > upper) {
        return 0.0;
    }
    return 1.0/(upper-lower);
}


// Metropolis-Hastings MCMC sampler of the UV luminosity fit likelihood
vector<vector<double> > mcmc_sampling(cosmology &C, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int Ns, int Nbi, rgen &mt) {
    vector<vector<double> > chain;
    vector<double> element;
        
    // Create normal distributions for each parameter based on its proposal width
    vector<normal_distribution<>> proposal_distributions;
    for (double step : steps) {
        proposal_distributions.push_back(normal_distribution<>(0.0, step));
    }
    
    vector<double> current = initial;
    double priorratio = 1.0;
    
    double logLcurrent = loglikelihood(C, data, current);
    cout << logLcurrent << endl;
    
    element.insert(element.end(), current.begin(), current.end());
    element.push_back(logLcurrent);
    chain.push_back(element);
    element.clear();
    
    double logLproposed, paccept, acceptanceratio = 0.0;
    vector<double> proposed;
    for (int i = 0; i < Ns; i++) {
        // generate new point by modifying the current sample
        proposed = current;
        for (int j = 0; j < proposed.size(); j++) {
            proposed[j] += proposal_distributions[j](mt);
        }
        
        priorratio = 1.0;
        for (int j = 0; j < priors.size(); ++j) {
            priorratio *= prior(proposed[j],priors[j])/prior(current[j],priors[j]);
        }
        
        if (priorratio > 0.0) {
            // calculate the acceptance ratio based on the likelihood
            logLproposed = loglikelihood(C, data, proposed);
            paccept = randomreal(0.0,1.0,mt);
                        
            // accept or reject the proposed sample
            if (logLproposed - logLcurrent - log(priorratio) > log(paccept)) {
                
                current = proposed;
                logLcurrent = logLproposed;
                cout << logLcurrent << endl;
                
                if (i > Nbi) {
                    acceptanceratio += 1.0/(1.0*(Ns-Nbi));
                    
                    element.insert(element.end(), current.begin(), current.end());
                    element.push_back(logLcurrent);
                    chain.push_back(element);
                    element.clear();
                }
            }
        }
    }
    cout << "acceptance ratio = " << acceptanceratio << endl;
    
    return chain;
}


vector<double> UFfit(cosmology &C, vector<vector<double> > &priors, vector<double> &steps, int Nsteps, int Nbi, int Nchains) {
    // random number generator
    rgen mt(time(NULL));
    
    // read UV luminosity data files, {MUV, Phi, +sigmaPhi, -sigmaPhi}
    vector<string> datafiles {"UVLF_2102.07775.txt", "UVLF_2108.01090.txt", "UVLF_2403.03171.txt"};
    vector<vector<double> > data = readUVdata(datafiles);
    
    // output the MCMC chains and find the best fit
    ofstream outfile;
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
    
    // return the best fit point
    return bf;
}



