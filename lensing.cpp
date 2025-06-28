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


// lensing convergence and its average, {kappa(r), bar{kappa}(<r)}
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


// MCMC sampling function
vector<vector<double> > MH2D(int N, int burnin, function<double(vector<double>&)> pdf, vector<double> &initial, vector<double> &steps, vector<vector<double> > &priors, function<double(vector<double>&)> cut, rgen &mt) {
    
    vector<vector<double> > samples(N);
    
    int Npar = initial.size();
    
    // create normal distributions for each parameter based on its proposal width
    vector<normal_distribution<>> proposal_distributions;
    for (double step : steps) {
        proposal_distributions.push_back(normal_distribution<>(0.0, step));
    }
    
    vector<double> current = initial;
    double pcurrent = pdf(current);
    
    vector<double> prop;
    double pprop, priorratio, paccept;
    for (int j = -burnin; j < N;) {
                
        // propose new point
        prop = current;
        for (int i = 0; i < Npar; i++) {
            prop[i] += proposal_distributions[i](mt);
        }
        
        // compute prior ratio (0 or 1, flat priors)
        priorratio = 1.0;
        for (int i = 0; i < Npar; i++) {
            priorratio *= prior(prop[i], priors[i]);
        }
        if (priorratio > 0.0) {
            priorratio = cut(prop);
        }
        
        // accept/reject
        if (priorratio > 0.0) {
            pprop = pdf(prop);
            paccept = randomreal(0.0,1.0,mt);
                                    
            if (paccept < min(1.0, pprop/pcurrent)) {
                current = prop;
                pcurrent = pprop;
                if (j >= 0) {
                    samples[j] = current;
                }
                j++;
            }
        }
    }
    
    // shuffle the samples
    shuffle(samples.begin(), samples.end(), mt);
    return samples;
}


// probability distribution of lnmu, {lnmu, dP/dlnmu}
vector<vector<double> > Plnmuf(cosmology &C, double zs, double rS, int Nhalos, int Nreal, int Nbins, rgen &mt) {
    
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
    
    double Nbar = Nhf(C, zs, rS, kappathr);
    poisson_distribution<int> PN(Nbar);
    
    cout << kappathr << "   " << Nbar << endl;
    
    // 2D PDF of z_l and ln M_l
    function<double(vector<double>&)> pdf = [&C, zs, rS, kappathr](vector<double> &zlnM) {
        double zl = zlnM[0];
        double M = exp(zlnM[1]);
        double rmax = rmaxf(C, zs, zl, M, rS, kappathr);
        return pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*interpolate2(zl, M, C.zlist, C.Mlist, C.HMFlist)[0];
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
    vector<double> gammalist(Nreal, 0.0);
    vector<double> lnmu;
    
    // generate the first list of z_l and M values
    int burnin = 1e3;
    int Nrand = 1e4;
    vector<double> initial = {zs/2.0, log(1.0e12)};
    vector<vector<double> > zlnMlist = MH2D(Nrand, burnin, pdf, initial, steps, priors, cut, mt);
    //writeToFile(zlnMlist, "zlnMlist.dat");
    
    int Nh, idx = 0;
    vector<double> kappav;
    double zl, M, rmax, r, phi, kappaj, gamma1j, gamma2j, muj = 1.0, lnmumax = 0.0, lnmumin = 0.0;
    double meankappa = 0.0, meangamma1 = 0.0, meangamma2 = 0.0;
    for (int j = 0; j < Nreal; j++) {
        kappaj = 0.0;
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
                zlnMlist = MH2D(Nrand, burnin, pdf, initial, steps, priors, cut, mt);
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
        meangamma1 += gamma1j;
        meangamma2 += gamma2j;
        
        kappalist[j] = kappaj;
        gamma1list[j] = gamma1j;
        gamma2list[j] = gamma2j;
    }
    meankappa = meankappa/(1.0*Nreal);
     
    double Alnmu = 0.0, Alnmu2 = 0.0;
    for (int j = 0; j < Nreal; j++) {
        kappalist[j] = kappalist[j] - meankappa;
        gammalist[j] = sqrt(pow(gamma1list[j], 2.0) + pow(gamma2list[j], 2.0));
        
        muj = 1.0/(pow(1.0-kappalist[j], 2.0) - pow(gammalist[j], 2.0));
        if (muj > 0.0) {
            lnmumin = min(log(muj), lnmumin);
            lnmumax = max(log(muj), lnmumax);
            Alnmu += log(muj);
            Alnmu2 += pow(log(muj),2.0);
            
            lnmu.push_back(log(muj));
        }
    }
    Alnmu = Alnmu/(1.0*Nreal);
    Alnmu2 = Alnmu2/(1.0*Nreal);
    double varlnmu = Alnmu2 - pow(Alnmu,2.0);
    
    // writeToFile(kappalist,"kappalist.dat");
    // writeToFile(gammalist,"gammalist.dat");
    
    // binning
    double dlnmu = sqrt(varlnmu)/(1.0*(Nbins-1));
    int Nbins2 = (int) ceil((lnmumax-lnmumin)/dlnmu + 1.0);
    vector<vector<double> > Plnmu(Nbins2, vector<double> (2, 0.0));
    
    for (int j = 0; j < Nbins2; j++) {
        Plnmu[j][0] = lnmumin + j*dlnmu;
    }
    int jb;
    for (int j = 0; j < lnmu.size(); j++) {
        jb = (int) round((Nbins2-1)*(lnmu[j]-lnmumin)/(lnmumax-lnmumin));
        if (jb >= 0 && jb < Nbins2) {
            Plnmu[jb][1] += 1.0;
        }
    }
    
    // convert from image plane to source plane (P_S ~ P_I/mu) and normalize
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
