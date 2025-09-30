#include "cosmology.h"
#include "lensing.h"


// NFW lens surface density at radius x
double Sigmaf(double rs, double rhos, double x) {
    double f = 1.0/3.0;
    if (x > 1) {
        f = (1 - 2*atan(sqrt((x-1)/(1+x)))/sqrt(pow(x,2)-1))/(pow(x,2)-1);
    }
    if (x < 1) {
        f = (1 - 2*atanh(sqrt((1-x)/(1+x)))/sqrt(1-pow(x,2)))/(pow(x,2)-1);
    }
    return 2*rs*rhos*f;
}

// average of the NFW lens surface density within radius x
double barSigmaf(double rs, double rhos, double x) {
    double g = 1 + log(1.0/2.0);
    if (x > 1) {
        g = 2*atan(sqrt((x-1)/(1+x)))/sqrt(pow(x,2)-1) + log(x/2);
    }
    if (x < 1) {
        g = 2*atanh(sqrt((1-x)/(1+x)))/sqrt(1-pow(x,2)) + log(x/2);
    }
    return 4*rs*rhos*g/pow(x,2);
}

// derivative of the NFW lens surface density at radius x
double DSigmaf(double rs, double rhos, double x) {
    double f = 0.0;
    if (x > 1) {
        f = -(1+2*pow(x,2))/(x*pow(pow(x,2)-1,2)) + 6*x*atan(sqrt((x-1)/(1+x)))/pow(pow(x,2)-1,5.0/2.0);
    }
    if (x < 1) {
        f = -(1+2*pow(x,2))/(x*pow(pow(x,2)-1,2)) + 6*x*atanh(sqrt((1-x)/(1+x)))/pow(1-pow(x,2),5.0/2.0);
    }
    return 2*rhos*f;
}

// lensing convergence, its average and its derivative, {kappa(r), bar{kappa}(<r), kappa'(r)}
// r_S = source size, r = impact parameter
vector<double> kappa(cosmology &C, double zs, double zl, double r, double M, double rS) {
    
    // angular diameter distances
    double DsA = C.DL(zs)/pow(1+zs,2.0);
    double DlA = C.DL(zl)/pow(1+zl,2.0);
    double DlsA = DsA - DlA*(1+zl)/(1+zs);
    
    double Sigmac = 2.08871e16*DsA/(4.0*PI*DlA*DlsA);
    
    // NFW scale radius and density
    vector<double> NFWparams = interpolate2(zl, M, C.zlist, C.Mlist, C.NFWlist);
    double rs = NFWparams[0];
    double rhos = NFWparams[1];
    
    double x, Sigma = 0.0, barSigma = 0.0, DSigma = 0.0;
    if (rS > 0.0) {
        // projection of the source to the lens plane
        double rSproj = rS*DlA/DsA;
        double Aproj = PI*pow(rSproj,2.0)/2.0;
        
        // average over the source projection
        int NR = 10, Ntheta = 10;
        double dR = rSproj/(1.0*(NR-1));
        double dtheta = PI/(1.0*(Ntheta-1));
        double R, theta;
        for (int jR = 0; jR < NR; jR++) {
            R = jR*dR;
            for (int jtheta = 0; jtheta < Ntheta; jtheta++) {
                theta = jtheta*dtheta;
                x = sqrt(pow(r,2.0) + pow(R,2.0) - 2.0*r*R*cos(theta))/rs;
                Sigma += Sigmaf(rs, rhos, x)*R*dR*dtheta;
                barSigma += barSigmaf(rs, rhos, x)*R*dR*dtheta;
                DSigma += DSigmaf(rs, rhos, x)*R*dR*dtheta;
            }
        }
        Sigma = Sigma/Aproj;
        barSigma = barSigma/Aproj;
        DSigma = DSigma/Aproj;
    }
    else {
        x = r/rs;
        Sigma = Sigmaf(rs, rhos, x);
        barSigma = barSigmaf(rs, rhos, x);
        DSigma = DSigmaf(rs, rhos, x);
    }
    return {Sigma/Sigmac, barSigma/Sigmac, DSigma/Sigmac};
}

// maximal r so that kappa^(1) > kappa_thr
double rmaxf(cosmology &C, double zs, double zl, double M, double rS, double kappathr) {
    double rmax;
    double logr1 = log(1.0e-6), logr2 = log(1.0e6);
    if (kappa(C, zs, zl, exp(logr1), M, rS)[0] > kappathr) {
        while (logr2-logr1 > 0.02) {
            rmax = exp((logr2+logr1)/2.0);
            if (kappa(C, zs, zl, rmax, M, rS)[0] > kappathr) {
                logr1 = log(rmax);
            } else {
                logr2 = log(rmax);
            }
        }
        rmax = exp((logr2+logr1)/2.0);
    } else {
        rmax = 0.0;
    }
    return rmax;
}

// number of halos with kappa^(1) > kappa_thr
double Nhf(cosmology &C, double zs, double rS, double kappathr) {
    double Nh = 0.0;
    double zl, dz, M, dlnM, dndlnM, rmax;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if(zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                rmax = rmaxf(C, zs, zl, M, rS, kappathr);
                Nh += 306.535*PI*pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
            }
        }
    }
    return Nh;
}

// variance of kappa from weak lenses
double sigmakappaW(cosmology &C, double zs, double rS, double kappathr) {
    double Nh = 0.0, kappa1 = 0.0, kappa2 = 0.0;
    double zl, dz, M, dlnM, dndlnM, r, kappar;
    int Nr = 100;
    double dlnr = 0.01;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if(zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                
                r = rmaxf(C, zs, zl, M, rS, kappathr);
                if (r == 0.0) {
                    r = 1.0e-6;
                }
                
                kappar = kappathr;
                while (kappar > 0.001*kappathr) {
                    kappar = kappa(C, zs, zl, r, M, rS)[0];
                    Nh += 306.535*PI*pow((1.0+zl)*r,2.0)/C.Hz(zl)*dndlnM*dlnr*dlnM*dz;
                    kappa1 += 306.535*PI*pow((1.0+zl)*r,2.0)/C.Hz(zl)*dndlnM*kappar*dlnr*dlnM*dz;
                    kappa2 += 306.535*PI*pow((1.0+zl)*r,2.0)/C.Hz(zl)*dndlnM*pow(kappar,2.0)*dlnr*dlnM*dz;
                    r = exp(log(r) + dlnr);
                }
            }
        }
    }
    return sqrt(kappa2 - pow(kappa1,2.0)/Nh);
}


// flat priors
double prior(double x, vector<double> &bounds) {
    double lower = bounds[0], upper = bounds[1];
    if (lower == upper) {
        return 1.0;
    }
    if (x < lower || x > upper) {
        return 0.0;
    }
    return 1.0/(upper-lower);
}


// sample N values from a PDF using Metropolis-Hastings MCMC sampler
vector<vector<double> > MCMC_sampling(int N, int Nburnin, function<double(vector<double>&)> logpdf, vector<double> &initial, vector<double> &steps, vector<vector<double> > &priors, function<double(vector<double>&)> cut, rgen &mt, int print, string filename) {
    
    ofstream outfile;
    if (print > 0) {
        outfile.open(filename);
    }
    
    vector<vector<double> > samples(N);
    int Npar = initial.size();
    
    // create normal distributions for each parameter based on its proposal width
    vector<normal_distribution<>> proposal_distributions;
    for (double step : steps) {
        proposal_distributions.push_back(normal_distribution<>(0.0, step));
    }
    
    vector<double> current = initial;
    double logpcurrent = logpdf(current);
    
    vector<double> prop;
    double logpprop, priorratio, paccept;
    int nnew = 0;
    for (int j = -Nburnin; j < N;) {
        
        // propose new point
        prop = current;
        for (int i = 0; i < Npar; i++) {
            prop[i] += proposal_distributions[i](mt);
        }
        
        // compute prior ratio
        priorratio = 1.0;
        for (int i = 0; i < Npar; i++) {
            priorratio *= prior(prop[i], priors[i])/prior(current[i], priors[i]);
        }
        
        // apply the cut
        if (priorratio > 0.0) {
            priorratio = cut(prop);
        }
        
        if (priorratio > 0.0) {
            logpprop = logpdf(prop);
            paccept = randomreal(0.0,1.0,mt);
            
            // accept or reject the proposed sample
            if (log(paccept) < logpprop - logpcurrent - log(priorratio)) {
                current = prop;
                logpcurrent = logpprop;
                if (j >= 0) {
                    nnew++;
                }
            }
            
            // after burn-in, add the element to the chain
            if (j >= 0) {
                if (print > 0) {
                    for (int jp = 0; jp < Npar; jp++) {
                        outfile << current[jp] << "   ";
                    }
                    outfile << endl;
                }
                samples[j] = current;
            }
            j++;
        }
        
        if (print > 0) {
            cout << j << "   " << nnew << "   " << "\r" << flush;
        }
    }
    
    if (print > 0) {
        outfile.close();
        cout << "acceptance ratio = " << nnew/(1.0*N) << endl;
    }
    
    return samples;
}


// probability distribution of lnmu, {lnmu, dP/dlnmu}
vector<vector<double> > Plnmuf(cosmology &C, double zs, double rS, int Nhalos, int Nreal, int Nbins, rgen &mt, int write) {
    
    // find threshold kappa that gives N_halos halos
    double kappa1 = 1.0e-12, kappa2 = 1.0;
    double kappathr = pow(10.0, (log10(kappa1) + log10(kappa2))/2.0);
    while (log10(kappa2) - log10(kappa1) > 0.01) {
        if (Nhf(C, zs, rS, kappathr) > Nhalos) {
            kappa1 = kappathr;
        } else {
            kappa2 = kappathr;
        }
        kappathr = pow(10.0, (log10(kappa1) + log10(kappa2))/2.0);
    }
    
    // distribution of kappa < kappa_thr
    double skappaW = sigmakappaW(C, zs, rS, kappathr);
    normal_distribution<double> PkappaW(0.0, skappaW);
    
    //cout << kappathr << "   " << skappaW << endl;
    if (skappaW < 0.0) {
        cout << "Error: negative standard deviation." << endl;
    }
    
    //double Nbar = Nhf(C, zs, rS, kappathr);
    double Nbar = 1.0*Nhalos;
    poisson_distribution<int> PN(Nbar);
    
    // 2D PDF of z_l and ln M_l
    function<double(vector<double>&)> logpdf = [&C, zs, rS, kappathr](vector<double> &zlnM) {
        double zl = zlnM[0];
        double M = exp(zlnM[1]);
        double rmax = rmaxf(C, zs, zl, M, rS, kappathr);
        return log(pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*interpolate2(zl, M, C.zlist, C.Mlist, C.HMFlist)[0]);
    };
    
    // prior in z_l - M_l plane so that kappa > kappa_thr
    function<double(vector<double>&)> cut = [&C, zs, rS, kappathr](vector<double> &zlnM) {
        double zl = zlnM[0];
        double M = exp(zlnM[1]);
        double kappa0 = kappa(C, zs, zl, 1.0e-6, M, rS)[0];
        if (kappa0 > kappathr) {
            return 1.0;
        }
        return 0.0;
    };
        
    // priors for z_l and ln M
    vector<vector<double> > priors = {{C.zmin, zs}, {log(C.Mmin), log(C.Mmax)}};
    vector<double> steps(priors.size());
    for (int j = 0; j < steps.size(); j++) {
        steps[j] = (priors[j][1] - priors[j][0])/10.0;
    }
    
    vector<double> kappalist(Nreal, 0.0);
    vector<double> gamma1list(Nreal, 0.0);
    vector<double> gamma2list(Nreal, 0.0);
    
    // generate the first list of z_l and M values
    int burnin = 1e3;
    int Nrand = 1e4;
    vector<double> initial = {zs/2.0, log(1.0e12)};
    vector<vector<double> > zlnMlist = MCMC_sampling(Nrand, burnin, logpdf, initial, steps, priors, cut, mt, 0, "");
    //writeToFile(zlnMlist, "zlnMlist.dat");
    
    int Nh, idx = 0;
    vector<double> kappav;
    double zl, M, rmax, r, phi, kappaj, gamma1j, gamma2j, gammaj, muj;
    double meankappa = 0.0;
    for (int j = 0; j < Nreal; j++) {
        kappaj =  PkappaW(mt); // include contribution of lenses with kappa < kappa_thr
        gamma1j = 0.0;
        gamma2j = 0.0;
        
        // pick number of halos from Poisson distribution
        Nh = PN(mt);
        
        // generate N_h halos
        for (int i = 0; i < Nh; i++) {
            zl = zlnMlist[idx][0];
            M = exp(zlnMlist[idx][1]);
            idx++;
                        
            if (idx >= Nrand) {
                // generate new set of z_l and M values
                initial = {zl, log(M)};
                zlnMlist = MCMC_sampling(Nrand, burnin, logpdf, initial, steps, priors, cut, mt, 0, "");
                idx = 0;
            }
            
            // find r_max so that kappa(r_max) = kappa_thr
            rmax = rmaxf(C, zs, zl, M, rS, kappathr);
                        
            r = sqrt(randomreal(0.0,1.0,mt))*rmax;
            phi = randomreal(0.0,2*PI,mt);
            
            kappav = kappa(C, zs, zl, r, M, rS);
            
            kappaj += kappav[0];
            gamma1j += -cos(2*phi)*(kappav[1]-kappav[0]);
            gamma2j += -sin(2*phi)*(kappav[1]-kappav[0]);
        }
        meankappa += kappaj;
        
        kappalist[j] = kappaj;
        gamma1list[j] = gamma1j;
        gamma2list[j] = gamma2j;
    }
    meankappa = meankappa/(1.0*Nreal);
    
    // compute mu
    vector<double> lnmulist;
    vector<double> lnmuAlist;
    vector<double> ln1pkappalist;
    vector<double> lngammalist;
    for (int j = 0; j < Nreal; j++) {
        kappaj = kappalist[j] - meankappa;
        gammaj = sqrt(pow(gamma1list[j], 2.0) + pow(gamma2list[j], 2.0));
        muj = 1.0/(pow(1.0-kappaj, 2.0) - pow(gammaj, 2.0));
        
        if (muj > 0.0 && gammaj > 0.0) {
            lnmulist.push_back(log(muj));
            lnmuAlist.push_back(log(1.0/pow(1.0-kappaj, 2.0)));
            ln1pkappalist.push_back(log(1.0+kappaj));
            lngammalist.push_back(log(gammaj));
        }
    }
    
    // binning
    vector<vector<double> > Plnmu = binSample(lnmulist, Nbins);
    if (write > 0) {
        vector<vector<double> > PlnmuA = binSample(lnmuAlist, Nbins);
        vector<vector<double> > Pln1pkappa = binSample(ln1pkappalist, Nbins);
        vector<vector<double> > Plngamma = binSample(lngammalist, Nbins);
            
        writeToFile(Plnmu,"Plnmu.dat");
        writeToFile(PlnmuA,"PlnmuA.dat");
        writeToFile(Pln1pkappa,"Pln1pkappa.dat");
        writeToFile(Plngamma,"Plngamma.dat");
    }
    
    // convert from image plane to source plane (P_S ~ P_I/mu) and normalize
    double dlnmu = Plnmu[1][0] - Plnmu[0][0];
    double norm = 0.0;
    for (int j = 0; j < Plnmu.size(); j++) {
        Plnmu[j][1] *= exp(-Plnmu[j][0]);
        norm += Plnmu[j][1]*dlnmu;
    }
    for (int j = 0; j < Plnmu.size(); j++) {
        Plnmu[j][1] *= 1.0/norm;
    }
    return Plnmu;
}


// loglikelihood of the Hubble digram data
double loglikelihood(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &par, int lensing, int dm, rgen &mt) {
    
    // initialize cosmology
    C.OmegaM = par[0];
    C.sigma8 = par[1];
    C.h = par[2];
    C.initialize(dm, pow(10.0, par[3]));
    
    double z, DL0, DL, sigmaDL, Y, dY, Pdet;
    
    // compute loglikelihood
    double logL = 0.0;
    if (lensing == 0) { // model without lensing
        for (int j = 0; j < data.size(); j++) {
            z = data[j][0];
            DL0 = C.DL(z);
            DL = data[j][1];
            sigmaDL = data[j][2];
            
            // compute P_det
            dY = sigmaDL/10.0;
            Y = min(DL0, DLthr) - 3.0*sigmaDL;
            Pdet = 0.0;
            while (Y <= min(DLthr, DL0 + 3.0*sigmaDL)) {
                Pdet += dY*NPDF(Y, DL0, sigmaDL);
                Y += dY;
            }
            logL += logNPDF(DL, DL0, sigmaDL) - log(Pdet);
        }
    }
    else { // model with lensing
        
        // compute the lensing distributions
        int Nreal = 2e3; // realizations
        int Nhalos = 40; // number of halos in each realization
        int Nbins = 12; // P(lnmu) bins
        double rS = 0.0; // lensing source radius in kpc
        
        vector<vector<vector<double> > > Plnmuz(C.Zlist.size());
        for (int jz = 0; jz < C.Zlist.size(); jz++) {
            z = C.Zlist[jz];
            Plnmuz[jz] = Plnmuf(C, z, rS, Nhalos, Nreal, Nbins, mt, 0);
        }
        
        int jz;
        vector<vector<double> > Plnmu;
        double L, dlnmu;
        for (int j = 0; j < data.size(); j++) {
            z = data[j][0];
            
            DL = data[j][1];
            sigmaDL = data[j][2];
            
            // find the closest z at which P(lnmu) is computed
            jz = lower_bound(C.Zlist.begin(), C.Zlist.end(), z)- C.Zlist.begin();
            if ((jz > 0 && C.Zlist[jz]-z > z-C.Zlist[jz-1]) || jz >= C.Zlist.size()) {
                jz--;
            }
            Plnmu = Plnmuz[jz];
            dlnmu = Plnmu[1][0] - Plnmu[0][0];
            
            // integrate over lnmu
            dY = sigmaDL/10.0;
            Pdet = 0.0, L = 0.0;
            for (int i = 0; i < Plnmu.size(); i++) {
                DL0 = C.DL(z)/exp(Plnmu[i][0]/2.0);
                
                L += dlnmu*Plnmu[i][1]*NPDF(DL, DL0, sigmaDL);
                
                // integrate over Y
                Y = min(DL0, DLthr) - 3.0*sigmaDL;
                while (Y <= min(DLthr, DL0 + 3.0*sigmaDL)) {
                    Pdet += dlnmu*Plnmu[i][1]*dY*NPDF(Y, DL0, sigmaDL);
                    Y += dY;
                }
            }
            
            logL += log(L/Pdet);
        }
    }
    
    return logL;
}


// MCMC inference of the Hubble diagram data
void Hubble_diagram_fit(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int N, int Nburnin, int lensing, int dm, rgen &mt, string filename) {
    
    // loglikelihood
    function<double(vector<double>&)> pdf = [&C, DLthr, &data, lensing, dm, &mt](vector<double> &par) {
        return loglikelihood(C, DLthr, data, par, lensing, dm, mt);
    };
    
    // no cut
    function<double(vector<double>&)> cut = [](vector<double> &par) {
        return 1.0;
    };
    
    MCMC_sampling(N, Nburnin, pdf, initial, steps, priors, cut, mt, 1, filename);
}
