#include "cosmology.h"
#include "lensing.h"


double Sigmacf(cosmology &C, double zs, double zl) {
    // angular diameter distances
    double DsA = C.DL(zs)/pow(1+zs,2.0);
    double DlA = C.DL(zl)/pow(1+zl,2.0);
    double DlsA = DsA - DlA*(1+zl)/(1+zs);
    
    // Sigma_c
    return 2.08871e16*DsA/(4.0*PI*DlA*DlsA);
}


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                              NFW halos                                                                         */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// from astro-ph/0608153
double FNFW(double x) {
    if (x > 1) {
        return (1 - 2*atan(sqrt((x-1)/(1+x)))/sqrt(pow(x,2)-1))/(pow(x,2)-1);
    }
    if (x < 1) {
        return (1 - 2*atanh(sqrt((1-x)/(1+x)))/sqrt(1-pow(x,2)))/(pow(x,2)-1);
    }
    return 1.0/3.0;
}
double gNFW(double x) {
    if (x > 1) {
        return 2*atan(sqrt((x-1)/(1+x)))/sqrt(pow(x,2)-1) + log(x/2);
    }
    if (x < 1) {
        return 2*atanh(sqrt((1-x)/(1+x)))/sqrt(1-pow(x,2)) + log(x/2);
    }
    return 1 + log(1.0/2.0);
}
double kappaNFW(double kappa0, double x) {
    return 2.0*kappa0*FNFW(x);
}
double gammaNFW(double kappa0, double x) {
    return 2.0*kappa0*(2.0*gNFW(x)/pow(x,2.0) - FNFW(x));
}
double kappa0NFW(double rs, double rhos, double Sigmac) {
    return rs*rhos/Sigmac;
}

// fix ellipticity using the fit of astro-ph/0508497
double epsilonNFW(cosmology &C, double z, double M) {
    double s = 0.54*pow(M/exp(interpolate(z, C.logMcharlist)), -0.05);
    return max(0.0, (1.0-s)/(1.0+s));
}

// kappa and gamma for pseudo elliptical NFW halos
vector<double> kappagammaNFWeps(double epsilon, double kappa0, double x, double phi) {
    double a1eps = 1.0-epsilon;
    double a2eps = 1.0+epsilon;
    double x1eps = sqrt(a1eps)*cos(phi)*x;
    double x2eps = sqrt(a2eps)*sin(phi)*x;
    double xeps = sqrt(pow(x1eps,2.0) + pow(x2eps,2.0));
    double phieps = atan(x2eps/x1eps);
    
    double kappaeps0 = kappaNFW(kappa0, xeps);
    double gammaeps0 = gammaNFW(kappa0, xeps);

    double kappaeps = kappaeps0 + epsilon*cos(2.0*phieps)*gammaeps0;
    double gammaeps = sqrt(pow(gammaeps0,2.0) + 2.0*epsilon*cos(2.0*phieps)*gammaeps0*kappaeps0 + pow(epsilon,2.0)*(pow(kappaeps0,2.0) - pow(cos(2.0*phieps)*gammaeps0,2.0)));
    
    return {kappaeps, gammaeps};
}
vector<double> kappagammaNFW(cosmology &C, double zs, double zl, double r, double M, double phi, double epsilon) {
    double Sigmac = Sigmacf(C, zs, zl);
    
    // NFW scale radius and density
    vector<double> NFWparams = interpolate2(zl, M, C.zlist, C.Mlist, C.NFWlist);
    double rs = NFWparams[0];
    double rhos = NFWparams[1];
    
    double kappa0 = kappa0NFW(rs, rhos, Sigmac);
    
    return kappagammaNFWeps(epsilon, kappa0, r/rs, phi);
}

// maximal r so that kappa_NFW > kappa_thr
double rmaxfNFW(cosmology &C, double zs, double zl, double M, double kappathr) {
    double rmax;
    double logr1 = log(1.0e-6), logr2 = log(1.0e6);
    if (kappagammaNFW(C, zs, zl, exp(logr1), M, 0.0, 0.0)[0] > kappathr) {
        while (logr2-logr1 > 0.02) {
            rmax = exp((logr2+logr1)/2.0);
            if (kappagammaNFW(C, zs, zl, rmax, M, 0.0, 0.0)[0] > kappathr) {
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

// number of halos with kappa_NFW > kappa_thr in each z and M bin
vector<vector<vector<double> > > deltaNhfNFW(cosmology &C, double zs, double kappathr) {
    vector<vector<vector<double> > > dNh(C.Nz, vector<vector<double> > (C.NM, vector<double> (3, 0.0)));
    double zl, dz, M, Mb, dlnM, dndlnM, rmax, sigma, sigmab;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                sigma = C.sigmalist[jM][1];
                
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                rmax = rmaxfNFW(C, zs, zl, M, kappathr);
                dNh[jz][jM][0] = 306.535*PI*pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
                
                Mb = 2.0*PI*pow(rmax,2.0)*(C.dc(zl)-C.dc(zl-dz))*C.rhoM0;
                sigmab = interpolate(Mb, C.sigmalist);
                
                // see Baumann (5.129)
                dNh[jz][jM][1] = C.Dg(zl)*sigmab*C.halobias(zl, sigma);
                
                dNh[jz][jM][2] = rmax;
            }
        }
    }
    return dNh;
}

// total number of halos with kappa_NFW > kappa_thr
double NhfNFW(cosmology &C, double zs, double kappathr) {
    double Nh = 0.0;
    double zl, dz, M, dlnM, dndlnM, rmax;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                rmax = rmaxfNFW(C, zs, zl, M, kappathr);
                Nh += 306.535*PI*pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
            }
        }
    }
    return Nh;
}

// variance of kappa from weak lenses
double sigmakappaW(cosmology &C, double zs, double kappathr) {
    double Nh = 0.0, kappa1 = 0.0, kappa2 = 0.0;
    double zl, dz, M, dlnM, dndlnM, r, kappar;
    int Nr = 100;
    double dlnr = 0.01;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                                
                r = rmaxfNFW(C, zs, zl, M, kappathr);
                if (r == 0.0) {
                    r = 1.0e-6;
                }
                
                kappar = kappathr;
                while (kappar > 0.001*kappathr) {
                    kappar = kappagammaNFW(C, zs, zl, r, M, 0.0, 0.0)[0];
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


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                           MCMC sampler                                                                         */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

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


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                    PDF of amplifications                                                                       */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// find threshold kappa
double findkappathr(int N, function<double(double)> Nf) {
    double kappa1 = 1.0e-12, kappa2 = 1.0;
    double kappathr = pow(10.0, (log10(kappa1) + log10(kappa2))/2.0);
    while (log10(kappa2) - log10(kappa1) > 0.01) {
        if (Nf(kappathr) > N) {
            kappa1 = kappathr;
        } else {
            kappa2 = kappathr;
        }
        kappathr = pow(10.0, (log10(kappa1) + log10(kappa2))/2.0);
    }
    return kappathr;
}


// probability distribution of lnmu, {lnmu, dP/dlnmu}
vector<vector<double> > Plnmuf(cosmology &C, double zs, int Nhalos, int Nreal, int Nbins, rgen &mt, int bias, int ell, int write) {
    
    // fix threshold kappa
    function<double(double)> NfNFW = [&C, zs](double kappa) {
        return NhfNFW(C, zs, kappa);
    };
    double kappathrH = findkappathr(Nhalos, NfNFW);
        
    // distribution of kappa_NFW < kappa_thr
    double skappaW = sigmakappaW(C, zs, kappathrH);
    normal_distribution<double> PkappaW(0.0, skappaW);
    if (skappaW < 0.0) {
        cout << "Error: negative standard deviation." << endl;
    }
    
    vector<double> kappalist(Nreal, 0.0);
    vector<double> gamma1list(Nreal, 0.0);
    vector<double> gamma2list(Nreal, 0.0);
    for (int j = 0; j < Nreal; j++) {
        kappalist[j] = PkappaW(mt);
    }
    
    vector<vector<vector<double> > > dNH = deltaNhfNFW(C, zs, kappathrH);
    if (write > 0) {
        writeToFile(C.zlist, C.Mlist, dNH, "dNH.dat");
    }
    
    vector<double> kappagamma;
    normal_distribution<double> pG(1.0);
    poisson_distribution<int> PN;
    double zl, M, rmax, r, phi, phiH, epsilon = 0.0, barN, sigma, deltab, lambda, meankappa = 0.0;
    int NH = 0;
    for (int jz = 0; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        if (zl < zs) {
            for (int jM = 0; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                
                // generate realizations of halos
                barN = dNH[jz][jM][0];
                sigma = dNH[jz][jM][1];
                rmax = dNH[jz][jM][2];
                for (int j = 0; j < Nreal; j++) {
                    
                    // bias
                    deltab = sigma*pG(mt);
                    lambda = barN*(1.0+deltab);
                    
                    if (bias == 0) {
                        lambda = barN;
                    }
                    
                    if (lambda < 0.2) { // if lambda is small, compare lambda to a random number U(0,1) (faster)
                        if (lambda > randomreal(0.0, 1.0, mt)) {
                            if (ell > 0) {
                                epsilon = epsilonNFW(C, zl, M); // pseudo ellipsity of the halo
                            }
                            r = sqrt(randomreal(0.0,1.0,mt))*rmax; // distance from the line-of-sight
                            phi = randomreal(0.0,2*PI,mt); // polar angle of r vector
                            phiH = randomreal(0.0,2*PI,mt); // orientation of the halo ellipticity
                            
                            kappagamma = kappagammaNFW(C, zs, zl, r, M, phiH, epsilon);
                            
                            kappalist[j] += kappagamma[0];
                            gamma1list[j] += cos(phi)*kappagamma[1];
                            gamma2list[j] += sin(phi)*kappagamma[1];
                            
                            meankappa += kappagamma[0];
                        }
                    }
                    else { // for larger lambda, generate number of halos from Poisson distribution (slower)
                        PN = poisson_distribution<int>(lambda);
                        NH = PN(mt);
                        if (NH > 0) {
                            if (ell > 0) {
                                epsilon = epsilonNFW(C, zl, M); // pseudo ellipsity of the halo
                            }
                            for (int jH = 0; jH < NH; jH++) {
                                r = sqrt(randomreal(0.0,1.0,mt))*rmax; // distance from the line-of-sight
                                phi = randomreal(0.0,2*PI,mt); // polar angle of r vector
                                phiH = randomreal(0.0,2*PI,mt); // orientation of the halo ellipticity
                                
                                kappagamma = kappagammaNFW(C, zs, zl, r, M, phiH, epsilon);
                                
                                kappalist[j] += kappagamma[0];
                                gamma1list[j] += cos(phi)*kappagamma[1];
                                gamma2list[j] += sin(phi)*kappagamma[1];
                                
                                meankappa += kappagamma[0];
                            }
                        }
                    }
                }
            }
        }
    }
    meankappa = meankappa/(1.0*Nreal);
    
    // compute mu
    vector<double> lnmulist, lnmuAlist, ln1pkappalist, lngammalist;
    double kappaj, gamma1j, gamma2j, gammaj, muj;
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
            
        writeToFile(Plnmu,"Plnmu_z=" + to_string_prec(zs,1) + ".dat");
        writeToFile(PlnmuA,"PlnmuA_z=" + to_string_prec(zs,1) + ".dat");
        writeToFile(Pln1pkappa,"Pln1pkappa_z=" + to_string_prec(zs,1) + ".dat");
        writeToFile(Plngamma,"Plngamma_z=" + to_string_prec(zs,1) + ".dat");
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


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                      Likelihood analysis of the Hubble diagram                                                                 */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

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
        int Nreal = 4e3; // realizations
        int Nhalos = 40; // number of halos in each realization
        int Nbins = 12; // P(lnmu) bins
        
        vector<vector<vector<double> > > Plnmuz(C.Zlist.size());
        for (int jz = 0; jz < C.Zlist.size(); jz++) {
            z = C.Zlist[jz];
            Plnmuz[jz] = Plnmuf(C, z, Nhalos, Nreal, Nbins, mt, 1, 1, 0);
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
