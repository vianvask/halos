#include "cosmology.h"

// generates a list in log scale
vector<double> cosmology::loglist(double xmin, double xmax, int Nx) {
    if (Nx > 1) {
        double dlogx = (log(xmax)-log(xmin))/(1.0*(Nx-1));
        double x = xmin;
        vector<double> tmp(Nx,0.0);
        for (int j = 0; j < Nx; j++) {
            tmp[j] = x;
            x = exp(log(x) + dlogx);
        }
        return tmp;
    } else {
        vector<double> tmp;
        return tmp;
    }
}


// CDM matter transfer function (astro-ph/9709112)
double cosmology::TM (double k) {
    double keq = 0.00326227*Hz(zeq)/(1.0+zeq);
    double ksilk = 0.0016*pow(OmegaB*pow(h,2.0),0.52)*pow(OmegaM*pow(h,2.0),0.73)*(1+pow(10.4*OmegaM*pow(h,2.0),-0.95));
    double a1 = pow(46.9*OmegaM*pow(h,2.0),0.67)*(1+pow(32.1*OmegaM*pow(h,2.0),-0.532));
    double a2 = pow(12.0*OmegaM*pow(h,2.0),0.424)*(1+pow(45.0*OmegaM*pow(h,2.0),-0.582));
    double b1 = 0.944/(1+pow(458*OmegaM*pow(h,2.0),-0.708));
    double b2 = pow(0.395*OmegaM*pow(h,2.0),-0.026);
    double acnum = pow(a1,-OmegaB/OmegaM)*pow(a2,-pow(OmegaB/OmegaM,3.0));
    double bcnum = 1.0/(1+b1*(pow(OmegaC/OmegaM,b2)-1));
    double s2 = 44.5*1000*log(9.83/(OmegaM*pow(h,2.0)))/sqrt(1+10*pow(OmegaB*pow(h,2.0),3.0/4.0));
    double b3 = 0.313*pow(OmegaM*pow(h,2.0),-0.419)*(1+0.607*pow(OmegaM*pow(h,2.0),0.674));
    double b4 = 0.238*pow(OmegaM*pow(h,2.0),0.223);
    double zd = 1291*pow(OmegaM*pow(h,2.0),0.251)/(1+0.659*pow(OmegaM*pow(h,2.0),0.828))*(1+b3*pow(OmegaB*pow(h,2.0),b4));
    double Rd = 31.5*OmegaB*pow(h,2.0)*pow(T0/2.7,-4.0)/(zd/1000);
    function<double(double)> g2 = [&](double y) {
        return y*(-6*sqrt(1+y) + (2.0+3.0*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    };
    double ab = 2.07*keq*s2*pow(1+Rd,-3.0/4.0)*g2((1+zeq)/(1+zd));
    double bb = 0.5 + OmegaB/OmegaM + (3.0-2.0*OmegaB/OmegaM)*sqrt(1+pow(17.2*OmegaM*pow(h,2.0),2.0));
    double bnode = 8.41*pow(OmegaM*pow(h,2.0),0.435);
    function<double(double)> qk = [&](double k) {
        return k/(13.41*keq);
    };
    function<double(double)> fk = [&](double k) {
        return 1.0/(1+pow(k*s2/5.4,4.0));
    };
    function<double(double, double)> C1 = [&](double k, double ac) {
        return 14.2/ac + 386.0/(1+69.9*pow(qk(k),1.08));
    };
    function<double(double, double, double)> To1 = [&](double k, double ac, double bc) {
        return log(exp(1.0)+1.8*bc*qk(k))/(log(exp(1.0)+1.8*bc*qk(k))+C1(k,ac)*pow(qk(k),2.0));;
    };
    function<double(double)> jo = [&](double x) {
        return sin(x)/x;
    };
    function<double(double)> s3 = [&](double k) {
        return s2/pow(1+pow(bnode/(k*s2),3.0),1.0/3.0);
    };
    function<double(double)> TC = [&](double k) {
        return fk(k)*To1(k,1.0,bcnum)+(1-fk(k))*To1(k,acnum,bcnum);
    };
    function<double(double)> TB = [&](double k) {
        return (To1(k,1.0,1.0)/(1+pow(k*s2/5.2,2.0)) + ab/(1+pow(bb/(k*s2),3.0))*exp(-pow(k/ksilk,1.4)))*jo(k*s3(k));
    };
    return OmegaB/OmegaM*TB(k) + OmegaC/OmegaM*TC(k);
}


// variance of the CDM matter fluctuations
vector<double> cosmology::sigmaC(double M, double deltaH) {
    double RM = pow(3.0*M/(4.0*PI*rhoM0),1.0/3.0);
    double DRM = RM/(3.0*M);
    
    double kmax = 1000.0/RM;
    double kmin = 1.0e-6*kmax;
    double dlogk = (log(kmax)-log(kmin))/(1.0*(Nk-1));
    
    double sigma2 = 0.0, dsigma2 = 0.0;
    double k1, k2 = kmin;
    for (int jk = 0; jk < Nk; jk++) {
        k1 = k2;
        k2 = exp(log(k2)+dlogk);
        sigma2 += (k2-k1)*(pow(Ws(k1*RM)*Deltak(k1,deltaH),2.0)/k1 + pow(Ws(k2*RM)*Deltak(k2,deltaH),2.0)/k2)/2.0;
        dsigma2 += (k2-k1)*(2.0*k2*DRM*DWs(k2*RM)*Ws(k2*RM)*pow(Deltak(k2,deltaH),2.0)/k2 + 2.0*k1*DRM*DWs(k1*RM)*Ws(k1*RM)*pow(Deltak(k1,deltaH),2.0)/k1)/2.0;
    }
    return {sqrt(sigma2), dsigma2/(2.0*sqrt(sigma2))};
}


// variance of the FDM matter fluctuations
vector<double> cosmology::sigmaF(double M, double deltaH, double m22) {
    double RM = pow(3.0*M/(4.0*PI*rhoM0),1.0/3.0);
    double DRM = RM/(3.0*M);
    
    double kmax = min(1000.0/RM, 10.0*km22(m22));
    double kmin = 1.0e-6*kmax;
    double dlogk = (log(kmax)-log(kmin))/(1.0*(Nk-1));
    
    double sigma2 = 0.0, dsigma2 = 0.0;
    double k1, k2 = kmin;
    for (int jk = 0; jk < Nk; jk++) {
        k1 = k2;
        k2 = exp(log(k2)+dlogk);
        sigma2 += (k2-k1)*(pow(Ws(k1*RM)*DeltakF(k1,deltaH,m22),2.0)/k1 + pow(Ws(k2*RM)*DeltakF(k2,deltaH,m22),2.0)/k2)/2.0;
        dsigma2 += (k2-k1)*(2.0*k2*DRM*DWs(k2*RM)*Ws(k2*RM)*pow(DeltakF(k2,deltaH,m22),2.0)/k2 + 2.0*k1*DRM*DWs(k1*RM)*Ws(k1*RM)*pow(DeltakF(k1,deltaH,m22),2.0)/k1)/2.0;
    }
    return {sqrt(sigma2), dsigma2/(2.0*sqrt(sigma2))};
}


// variance of the WDM matter fluctuations
vector<double> cosmology::sigmaW(double M, double deltaH, double m3) {
    double RM = pow(3.0*M/(4.0*PI*rhoM0),1.0/3.0);
    double DRM = RM/(3.0*M);
    
    double kmax = min(1000.0/RM, 10.0*km3(m3));
    double kmin = 1.0e-6*kmax;
    double dlogk = (log(kmax)-log(kmin))/(1.0*(Nk-1));
    
    double sigma2 = 0.0, dsigma2 = 0.0;
    double k1, k2 = kmin;
    for (int jk = 0; jk < Nk; jk++) {
        k1 = k2;
        k2 = exp(log(k2)+dlogk);
        sigma2 += (k2-k1)*(pow(Ws(k1*RM)*DeltakW(k1,deltaH,m3),2.0)/k1 + pow(Ws(k2*RM)*DeltakW(k2,deltaH,m3),2.0)/k2)/2.0;
        dsigma2 += (k2-k1)*(2.0*k2*DRM*DWs(k2*RM)*Ws(k2*RM)*pow(DeltakW(k2,deltaH,m3),2.0)/k2 + 2.0*k1*DRM*DWs(k1*RM)*Ws(k1*RM)*pow(DeltakW(k1,deltaH,m3),2.0)/k1)/2.0;
    }
    return {sqrt(sigma2), dsigma2/(2.0*sqrt(sigma2))};
}

// variance of the enhanced matter fluctuations
vector<double> cosmology::sigmaE(double M, double deltaH, double kc) {
    double RM = pow(3.0*M/(4.0*PI*rhoM0),1.0/3.0);
    double DRM = RM/(3.0*M);
    
    double kmax = 1000.0/RM;
    double kmin = 1.0e-6*kmax;
    double dlogk = (log(kmax)-log(kmin))/(1.0*(Nk-1));
    
    double sigma2 = 0.0, dsigma2 = 0.0;
    double k1, k2 = kmin;
    for (int jk = 0; jk < Nk; jk++) {
        k1 = k2;
        k2 = exp(log(k2)+dlogk);
        sigma2 += (k2-k1)*(pow(Ws(k1*RM)*DeltakE(k1,deltaH,kc),2.0)/k1 + pow(Ws(k2*RM)*DeltakE(k2,deltaH,kc),2.0)/k2)/2.0;
        dsigma2 += (k2-k1)*(2.0*k2*DRM*DWs(k2*RM)*Ws(k2*RM)*pow(DeltakE(k2,deltaH,kc),2.0)/k2 + 2.0*k1*DRM*DWs(k1*RM)*Ws(k1*RM)*pow(DeltakE(k1,deltaH,kc),2.0)/k1)/2.0;
    }
    return {sqrt(sigma2), dsigma2/(2.0*sqrt(sigma2))};
}


// list of matter fluctuation variances, {M, sigma, dsigma/dM}
vector<vector<double> > cosmology::sigmalistf(double m22, double m3, double kc) {
    int Nextra = (int) ceil((NM-1)*log(3.0)/log(Mmax/Mmin)); // extend the M range for computation of dotM
    vector<vector<double> > Ms(NM+Nextra, vector<double> (3,0.0));
    
    double dlogM = (log(Mmax)-log(Mmin))/(1.0*(NM-1));
    double M = exp(log(Mmin) - dlogM);
        
    vector<double> sigma;
    for (int jM = 0; jM < NM+Nextra; jM++) {
        M = exp(log(M) + dlogM);
        if (m22 == 0.0 && m3 == 0.0) {
            sigma = sigmaC(M, deltaH8);
        }
        if (m22 > 0.0) {
            sigma = sigmaF(M, deltaH8, m22);
        }
        if (m3 > 0.0) {
            sigma = sigmaW(M, deltaH8, m3);
        }
        if (kc > 0.0) {
            sigma = sigmaE(M, deltaH8, kc);
        }
        Ms[jM][0] = M;
        Ms[jM][1] = sigma[0];
        Ms[jM][2] = sigma[1];
    }
    return Ms;
}

// first crossing probability with ellipsoidal collapse
double cosmology::pFC(double delta, double S) {
    double p = 0.3;
    double q = 0.8;
    double A = 1.0/(1+pow(2.0,-p)*tgammaf(0.5-p)/sqrt(PI));
    
    double nu2 = 0.0;
    if (S > 0.0) {
        nu2 = pow(delta,2.0)/S;
    }
    if (nu2 > 0.0) {
        return A*(1+pow(q*nu2,-p))*sqrt(q*nu2/(2.0*PI))*exp(-q*nu2/2.0)/S;
    }
    return 0.0;
}

// halo mass function and growth rate, {z, M, dn/dlnM, dotM, dotM/dM}
vector<vector<vector<double> > > cosmology::HMFlistf() {
    
    double dz = 0.01; // z step over which DeltaM is computed
    double q = 0.75; // for computation of dotM
    
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz-1));
    vector<vector<vector<double> > > dndlnM(Nz, vector<vector<double> > (NM, vector<double> (5,0.0)));
    
    double z, M, Mp1, S, S2, Ddeltac, dMdt, dMdtp1;
    vector<double> sigma, sigma2;
    for (int jz = 0; jz < Nz; jz++) {
        z = zlist[jz];
        
        sigma = sigmalist[0];
        S = pow(sigma[1],2.0);
        Ddeltac = (deltaell(z+dz,S)-deltaell(z,S))/dz;
        M = sigma[0];
        sigma2 = interpolaten(2.0*M, sigmalist);
        S2 = pow(sigma2[1],2.0);
        dMdt = (1+z)*Hz(z)*sqrt(2.0/PI)*Ddeltac/abs(2.0*sigma[1]*sigma[2])*sqrt(S-S2);
        
        for (int jM = 0; jM < NM; jM++) {
            dndlnM[jz][jM][0] = z;
            dndlnM[jz][jM][1] = M;
            dndlnM[jz][jM][2] = -rhoM0*pFC(deltac(z),pow(sigma[1],2.0))*2.0*sigma[1]*sigma[2];
            
            sigma = sigmalist[jM+1];
            S = pow(sigma[1],2.0);
            
            Ddeltac = (deltaell(z+dz,S)-deltaell(z,S))/dz;
            
            Mp1 = sigma[0];
            sigma2 = interpolaten(2.0*Mp1, sigmalist);
            S2 = pow(sigma2[1],2.0);
            dMdtp1 = (1+z)*Hz(z)*sqrt(2.0/PI)*Ddeltac/abs(2.0*sigma[1]*sigma[2])*sqrt(S-S2);
            
            dndlnM[jz][jM][3] = dMdt;
            dndlnM[jz][jM][4] = (dMdtp1 - dMdt)/(Mp1 - M);
            
            M = Mp1;
            dMdt = dMdtp1;
        }
    }
    return dndlnM;
}

// star formation rate f_*(M) and its derivative df_*/dM
double cosmology::fstar(double z, double M, double Mc, double Mt, double epsilon, double alpha, double beta) {
    if (alpha>0.0 && beta>0.0) {
        return epsilon*(alpha+beta)/(beta*pow(M/Mc,-alpha) + alpha*pow(M/Mc,beta))*exp(-Mt/M);
    }
    return epsilon;
}
double cosmology::Dfstarperfstar(double M, double Mc, double Mt, double epsilon, double alpha, double beta) {
    if (alpha>0.0 && beta>0.0) {
        return beta*((alpha+beta)/(alpha*pow(M/Mc,alpha+beta)+beta) - 1.0)/M + Mt/pow(M,2.0);
    }
    return 0.0;
}

// TODO: growth of stellar mass and BH mass by mergers
vector<double> cosmology::growbymergers(vector<double> &MJ, double z, double zp) {
    double M, Mp, dMp;
    vector<double> MJnew(NM, 0.0);
    vector<double> sigma, sigmap;
    /*for (int jM = NM-1; jM > 0; jM--) {
        sigma = sigmalist[jM];
        M = sigma[0];
        for (int jMp = 1; jMp < jM; jMp++) {
            sigmap = sigmalist[jMp];
            Mp = sigmap[0];
            dMp = Mp - sigmalist[jMp-1][0];
            
            MJnew[jM] += 2.0*sigmap[1]*abs(sigmap[2])*MJ[jMp]*M/Mp*pFC(deltac(zp),pow(sigmap[1],2.0),deltac(z),pow(sigma[1],2.0))*dMp;
        }
    }*/
    return MJnew;
}

// growth of stellar mass by star formation
vector<vector<double> > cosmology::evolvestellarmass(double Mc, double Mt, double epsilon, double alpha, double beta) {
    vector<vector<double> > Mst(Nz, vector<double> (NM, 0.0));
    
    double dt, M, dotM;
    double z, zp = zlist[Nz-1];
    
    for (int jz = Nz-2; jz >= 0; jz--) {
        z = zlist[jz];
        
        // star formation
        for (int jM = 0; jM < NM; jM++) {
            M = HMFlist[jz][jM][1];
            dotM = HMFlist[jz][jM][3];
            Mst[jz][jM] = Mst[jz+1][jM] + fstar(z,M,Mc,Mt,epsilon,alpha,beta)*dotM*(zp-z)/((1.0+z)*Hz(z));
        }
        
        // mergers
        Mst[jz] = growbymergers(Mst[jz], z, zp);
        
        zp = z;
    }
    return Mst;
}

// growth of the BH mass by accretion
vector<vector<double> > cosmology::evolveBHmass(double Mc, double Mt, double epsilon, double alpha, double beta, double fEdd, double facc1, double facc2) {
    vector<vector<double> > MBH(Nz, vector<double> (NM, 0.0));
    
    double dt, M, dotM, DeltaMEdd, frem;
    double z, dz, zp = zlist[Nz-1];

    for (int jz = Nz-2; jz >= 0; jz--) {
        z = zlist[jz];
        dz = zp - z;
        dt = dz/((1.0+z)*Hz(z));
        
        // accretion
        for (int jM = 0; jM < NM; jM++) {
            M = HMFlist[jz][jM][1];
            dotM = HMFlist[jz][jM][3];
            DeltaMEdd = 0.0022*MBH[jz][jM]*dt;
            if (M < Mc) {
                frem = 0.0;
            } else {
                frem = 1.0 - fstar(z,M,Mc,Mt,epsilon,alpha,beta)/epsilon;
            }
            MBH[jz][jM] = MBH[jz+1][jM] + min(frem*(facc1*dotM*dt + facc2*dt*M), fEdd*DeltaMEdd);
        }
        
        // mergers
        MBH[jz] = growbymergers(MBH[jz], z, dz);
        zp = z;
    }
    return MBH;
}


// halo consentration parameter (1601.02624)
double cosmology::cons(double sigma, double z) {
    double c0 = 3.395*pow(1+z,-0.215);
    double beta = 0.307*pow(1+z,0.540);
    double gamma1 = 0.628*pow(1+z,-0.047);
    double gamma2 = 0.317*pow(1+z,-0.893);
    double nu0 = 4.135 - 0.564*(1+z) - 0.210*pow(1+z,2.0) + 0.0557*pow(1+z,3.0) - 0.00348*pow(1+z,4.0);
    
    double nu = deltac(z)/sigma;
    
    return c0*pow(nu/nu0,-gamma1)*pow(1+pow(nu/nu0,1.0/beta),-beta*(gamma2-gamma1));
}

// halo consentration parameter, {z, M, c, dc/dM}
vector<vector<vector<double> > > cosmology::conslistf() {
    vector<vector<vector<double> > > zMc(Nz, vector<vector<double> > (NM, vector<double> (4,0.0)));
    double z, c, cn, M, Mn;
    for (int jz = 0; jz < Nz; jz++) {
        z = zlist[jz];
        M = sigmalist[0][0];
        c = cons(sigmalist[0][1],z);
        for (int jM = 0; jM < NM; jM++) {
            Mn = sigmalist[jM+1][0];
            cn = cons(sigmalist[jM+1][1],z);
            
            zMc[jz][jM][0] = z;
            zMc[jz][jM][1] = M;
            zMc[jz][jM][2] = c;
            zMc[jz][jM][3] = (cn - c)/(Mn - M);
            
            M = Mn;
            c = cn;
        }
    }
    return zMc;
}


// NFW scale radius and density and their derivatives, {r_s, dr_s/dM, rho_s, drho_s/dM}
vector<vector<vector<double> > > cosmology::NFWlistf() {
    vector<vector<vector<double> > > zMNFW(Nz, vector<vector<double> > (NM, vector<double> (6,0.0)));
    
    double z, M, c, Dc, r200, Dr200, rs, Drs, rhos, Drhos;
    vector<double> zMc(4,0.0);
    for (int jz = 0; jz < Nz; jz++) {
        for (int jM = 0; jM < NM; jM++) {
            
            zMc = conslist[jz][jM];
            z = zMc[0];
            M = zMc[1];
            c = zMc[2];
            Dc = zMc[3];
            
            r200 = pow(3.0*M/(4.0*PI*200*rhoc),1.0/3.0);
            Dr200 = 3.0/(4.0*PI*200*rhoc)*pow(3.0*M/(4.0*PI*200*rhoc),-2.0/3.0);
            
            rs =  r200/c;
            Drs =  Dr200/c - r200*Dc/pow(c,2.0);
            
            rhos = 200*rhoc*pow(c,3.0)*(1+c)/(3.0*((1+c)*log(1+c) - c));
            Drhos = Dc*pow(c,2.0)*(-c*(3+4*c) + 3*pow(1+c,2.0)*log(1+c))/(3.0*pow(c-(1+c)*log(1+c),2.0));
            
            zMNFW[jz][jM][0] = z;
            zMNFW[jz][jM][1] = M;
            zMNFW[jz][jM][2] = rs;
            zMNFW[jz][jM][3] = Drs;
            zMNFW[jz][jM][4] = rhos;
            zMNFW[jz][jM][5] = Drhos;
        }
    }
    return zMNFW;
}



// comoving distance, {z,d_c}
vector<vector<double> > cosmology::dclist() {
    int Nz2 = 100*Nz;
    
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz2-2));
    double z1 = 0.0, z2 = zmin;
    vector<double> d0 = {0.0, 0.0};
    
    vector<vector<double> > dlist(Nz2, vector<double> (2));
    dlist[0] = d0;
    
    for (int jz = 1; jz < Nz2; jz++) {
        d0[0] = z2;
        d0[1] += (z2-z1)*306.535*exp((log(1.0/Hz(z2))+log(1.0/Hz(z1)))/2.0);
        dlist[jz] = d0;
        
        z1 = z2;
        z2 = exp(log(z2) + dlogz);
    }
    return dlist;
}


// age of the Universe, {z,t}
vector<vector<double> > cosmology::tlist() {
    int Nz2 = 100*Nz;
    
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz2-2));
    double z1 = 0.0, z2 = zmin;
    vector<double> d0 = {0.0, 0.0};
    
    vector<vector<double> > tlist(Nz2, vector<double> (2));
    tlist[0] = d0;
    
    for (int jz = 1; jz < Nz2; jz++) {
        d0[0] = z2;
        d0[1] += (z2-z1)*exp((log(1.0/((1+z2)*Hz(z2)))+log(1.0/((1+z1)*Hz(z1))))/2.0);
        
        tlist[jz] = d0;
        
        z1 = z2;
        z2 = exp(log(z2)+dlogz);
    }
    
    double tmax = tlist[Nz2-1][1];
    for (int jz = 0; jz < Nz2; jz++) {
        tlist[jz][1] = tmax - tlist[jz][1];
    }
    
    return tlist;
}

