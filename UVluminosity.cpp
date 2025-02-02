#include "cosmology.h"
#include "UVluminosity.h"

double kappaUV = 1.15e-22; // Msun s erg^-1 Myr^-1

// star formation rate
double fstar(double M, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return epsilon*pow(pow(M/Mc,alpha/gamma) + pow(M/Mc,beta/gamma),-gamma);
}

// derivative of the star formation rate, df_*/dM
double Dfstar(double M, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return -(alpha/(1+pow(M/Mc,(-alpha+beta)/gamma)) + beta/(1+pow(M/Mc,(alpha-beta)/gamma)))*fstar(M,Mc,epsilon,alpha,beta,gamma)/M;
}

// the UV magnitude
double MUV(double M, double dotM, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return 51.63 - 1.08574*log(fstar(M,Mc,epsilon,alpha,beta,gamma)*dotM/kappaUV);
}
    
// derivative of the UV magnitude, dM_UV/dM
double DMUV(double M, double dotM, double DdotM, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return -1.08574*(DdotM/dotM + Dfstar(M,Mc,epsilon,alpha,beta,gamma)/fstar(M,Mc,epsilon,alpha,beta,gamma));
}

// dust extinction, MUV -> MUV - AUV (see 1406.1503 and Table 3 of 1306.2950)
double AUV(double MUV, double z) {
    double C0 = 4.4, C1 = 2.0, sigmabeta = 0.34;
    double beta = -exp(0.17*(19.5+MUV)/(1.54+0.075*z))*(1.54+0.075*z);
    return max(0.0, C0 + 0.2*log(10)*pow(C1*sigmabeta,2.0) + C1*beta);
}

// UV luminosity function as a function of halo mass M (see e.g. 1906.06296)
double UVLF(double M, double dotM, double DdotM, double dndlnM, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return -exp(-Mt/M)*dndlnM/M/DMUV(M,dotM,DdotM,Mc,epsilon,alpha,beta,gamma);
}

// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
vector<vector<vector<double> > > PhiUV(cosmology &C, vector<double> &Zlist, vector<vector<vector<double> > > &Plnmuz, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma) {
    
    vector<vector<vector<double> > > PhiUVlensed(Zlist.size(), vector<vector<double> > (C.NM, vector<double> (4,0.0)));
    vector<vector<double> > PhiUVlist(C.NM, vector<double> (3,0.0));
    
    int jzp;
    double z, dlnmu, Mj, dotMj, DdotMj, dndlnMj, MUVj, AUVj, PhiUVj;
    vector<vector<double> > Plnmu(Plnmuz[0].size(), vector<double> (2,0.0));
    for (int jZ = 0; jZ < Zlist.size(); jZ++) {
        z = Zlist[jZ];
        Plnmu = Plnmuz[jZ];
        
        // find index of z in the longer list of z values
        jzp = lower_bound(C.zlist.begin(), C.zlist.end(), z)- C.zlist.begin();
        if (jzp > 0 && C.zlist[jzp]-z > z-C.zlist[jzp-1]) {
            jzp--;
        }
             
        // compute the unlensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            Mj = C.dotMlist[jzp][jM][1];
            dotMj = C.dotMlist[jzp][jM][2];
            DdotMj = C.dotMlist[jzp][jM][3];
            dndlnMj = C.HMFlist[jzp][jM][2];
            
            MUVj = MUV(Mj, dotMj, Mc, epsilon, alpha, beta, gamma);
            
            PhiUVlist[jM][0] = MUVj;
            PhiUVlist[jM][1] = AUV(MUVj, z);
            PhiUVlist[jM][2] = UVLF(Mj, dotMj, DdotMj, dndlnMj, Mt, Mc, epsilon, alpha, beta, gamma);
        }
        
        // compute the lensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            MUVj = PhiUVlist[jM][0];
            AUVj = PhiUVlist[jM][1];
            
            // convolution integral
            PhiUVj = 0.0;
            dlnmu = Plnmu[1][0]-Plnmu[0][0];
            for (int jb = 0; jb < Plnmu.size(); jb++) {
                PhiUVj += interpolaten(MUVj - AUVj + 1.08574*Plnmu[jb][0], PhiUVlist)[2]*Plnmu[jb][1]*dlnmu/exp(Plnmu[jb][0]);
            }
            
            PhiUVlensed[jZ][jM][0] = MUVj;
            PhiUVlensed[jZ][jM][1] = PhiUVlist[jM][2]; // no dust + no lensing
            PhiUVlensed[jZ][jM][2] = interpolaten(MUVj - AUVj, PhiUVlist)[2]; // dust + lensing
            PhiUVlensed[jZ][jM][3] = PhiUVj; // dust + lensing
        }
    }
    return PhiUVlensed;
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
double NPDF2(double x, double mu, double sigmam, double sigmap) {
    if (x < mu) {
        return 2*sigmam*NPDF(x, mu, sigmam)/(sigmam+sigmap);
    }
    return 2*sigmap*NPDF(x, mu, sigmap)/(sigmam+sigmap);
}

// likelihood function
double loglikelihood(cosmology &C, vector<double> &Zlist, vector<vector<vector<double> > > &Plnmuz, vector<vector<double> > &data, vector<double> &params) {
    
    // compute the UV luminosity function
    vector<vector<vector<double> > > PhiUVlist = PhiUV(C, Zlist, Plnmuz, 1.0e9, params[0], params[1], params[2], params[3], 1.0);
    
    vector<vector<double> > PhiUVz;
    double z, MUV, meanPhi, sigmaPhip, sigmaPhim, Phi, logL = 0.0;
    for (int j = 0; j < data.size(); j++) {
        z = data[j][0];
        MUV = data[j][1];
        meanPhi = data[j][2];
        sigmaPhip = data[j][3];
        sigmaPhim = data[j][4];
        
        // z in Zlist
        PhiUVz = PhiUVlist[lower_bound(Zlist.begin(), Zlist.end(), z) - Zlist.begin()];
        
        // Phi(M_UV), convert to Mpc^-3 as data is in such units
        Phi = 1.0e9*interpolaten(MUV, PhiUVz)[3];
        
        logL += log(NPDF2(Phi, meanPhi, sigmaPhim, sigmaPhip));
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

// Metropolis-Hastings MCMC sampler
vector<vector<double> > mcmc_sampling(cosmology &C, vector<double> &Zlist, vector<vector<vector<double> > > &Plnmuz, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int Ns, rgen &mt) {
    vector<vector<double> > chain;
    vector<double> element;
        
    // Create normal distributions for each parameter based on its proposal width
    vector<normal_distribution<>> proposal_distributions;
    for (double step : steps) {
        proposal_distributions.push_back(normal_distribution<>(0.0, step));
    }
    
    vector<double> current = initial;
    double priorratio = 1.0;
    
    double logLcurrent = loglikelihood(C, Zlist, Plnmuz, data, current);
    cout << logLcurrent << endl;
    
    element.insert(element.end(), current.begin(), current.end());
    element.push_back(logLcurrent);
    chain.push_back(element);
    element.clear();
    
    double logLproposed;
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
            logLproposed = loglikelihood(C, Zlist, Plnmuz, data, proposed);
            
            // accept or reject the proposed sample
            if (logLproposed - logLcurrent - log(priorratio) > log(randomreal(0.0,1.0,mt))) {
                current = proposed;
                logLcurrent = logLproposed;
                cout << logLcurrent << endl;
                
                element.insert(element.end(), current.begin(), current.end());
                element.push_back(logLcurrent);
                chain.push_back(element);
                element.clear();
            }
        }
    }
    return chain;
}



