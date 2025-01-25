#include "cosmology.h"

// matter transfer function (astro-ph/9709112)
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


// variance of the matter fluctuations
double cosmology::sigmaf(double M, double deltaH) {
    double RM = pow(3.0*M/(4.0*PI*rhoM0),1.0/3.0);
    double kmin = 0.001/RM;
    double kmax = 1000.0/RM;
    double dlogk = (log(kmax)-log(kmin))/(1.0*(Nk-1));
    
    double sigma2 = 0.0;
    double k1, k2 = kmin;
    for (int jk = 0; jk < Nk; jk++) {
        k1 = k2;
        k2 = exp(log(k2)+dlogk);
        sigma2 += exp((log(pow(W(k2*RM)*Deltak(k2,deltaH),2.0)/k2) + log(pow(W(k1*RM)*Deltak(k1,deltaH),2.0)/k1))/2.0)*(k2-k1);
    }
    
    return sqrt(sigma2);
}
// variance of matter fluctuations, {M, sigma, dsigma/dM}
vector<vector<double> > cosmology::sigmalistf() {
    
    int Nextra = (int) ceil((NM-1)*log(3.0)/log(Mmax/Mmin));
        
    double dlogM = (log(Mmax)-log(Mmin))/(1.0*(NM-1));
    vector<vector<double> > Ms(NM+Nextra, vector<double> (3,0.0));
    double M = exp(log(Mmin) - dlogM);
    double sigma = sigmaf(M, deltaH8);
    double sigman, Mn;
    for (int jM = 0; jM < NM+Nextra; jM++) {
        Mn = exp(log(M) + dlogM);
        sigman = sigmaf(Mn, deltaH8);
        
        Ms[jM][0] = Mn;
        Ms[jM][1] = sigman;
        Ms[jM][2] = (sigman-sigma)/(Mn-M);

        M = Mn;
        sigma = sigman;
    }
    return Ms;
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
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz-1));
    vector<vector<vector<double> > > zMc(Nz, vector<vector<double> > (NM, vector<double> (4,0.0)));
    double z = zmin;
    double c, cn, M, Mn;
    for (int jz = 0; jz < Nz; jz++) {
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
        z = exp(log(z) + dlogz);
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


// Seth-Tormen HMF, {z,M,dn/dlnM}
vector<vector<vector<double> > > cosmology::hmflistf() {
    function<double(double)> nuf = [&](double nu) {
        double p = 0.3;
        double q = 0.75;
        double A = 1.0/(1+(pow(2.0,-p)*tgammaf(0.5-p)/sqrt(PI)));
        return A*(1+pow(q*nu,-p))*sqrt(q*nu/(2.0*PI))*exp(-q*nu/2.0);
    };
        
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz-1));
    vector<vector<vector<double> > > dndlnM(Nz, vector<vector<double> > (NM, vector<double> (3,0.0)));
    double z = zmin, M;
    for (int jz = 0; jz < Nz; jz++) {
        for (int jM = 0; jM < NM; jM++) {
            M = sigmalist[jM][0];
            
            dndlnM[jz][jM][0] = z;
            dndlnM[jz][jM][1] = M;
            dndlnM[jz][jM][2] = -rhoM0*nuf(pow(deltac(z)/sigmalist[jM][1],2.0))*2.0*sigmalist[jM][2]/sigmalist[jM][1];
        }
        z = exp(log(z) + dlogz);
    }
    
    return dndlnM;
}


// the probability that the halo whose mass at z' is M' ends up being a part of a halo in the mass range (M,M+dM) at z<z':
double cosmology::dPpdM(vector<double> &sigma, double z, vector<double> &sigmap, double zp) {
    return 2.0*sigma[1]*abs(sigma[2])*pow(sigmap[1]/sigma[1],3.0)*(deltac(zp)-deltac(z))*deltac(z)/deltac(zp)/sqrt(2.0*PI*pow(pow(sigmap[1],2.0)-pow(sigma[1],2.0),3.0))*exp(-pow(pow(sigmap[1],2.0)*deltac(z) - pow(sigma[1],2.0)*deltac(zp),2.0)/(2.0*pow(sigma[1]*sigmap[1],2.0)*(pow(sigmap[1],2.0)-pow(sigma[1],2.0))));
}
double cosmology::DeltaM(vector<double> &sigmap, double z, double dz) {
    double Mp = sigmap[0];
    
    // list of ratios at which the integrands are evaluated, M = (1+f)*M'
    vector<double> fM(71, 0.0);
    double logdf = 0.1;
    for (int jf = 0; jf < fM.size(); jf++) {
        fM[fM.size() - 1 - jf] = pow(10.0, -jf*logdf);
    }
    
    // integrals over M
    vector<double> sigmaj;
    double M, dM, dPdMj, dPdMjm1 = 0.0, num = 0.0, den = 0.0;
    for (int jM = 1; jM < fM.size(); jM++) {
        M = (1.0+fM[jM])*Mp;
        sigmaj = interpolaten(M, sigmalist);
        dPdMj = dPpdM(sigmaj, z, sigmap, z+dz);
        
        dM = M - (1.0+fM[jM-1])*Mp;
        num += M*(dPdMj + dPdMjm1)*dM/2.0;
        den += (dPdMj + dPdMjm1)*dM/2.0;
        
        dPdMjm1 = dPdMj;
    }
    return num/den - Mp;
}


// growth rate of the halo through mergers with smaller halos, {z,M,dM/dt}
vector<vector<vector<double> > > cosmology::dotMlistf() {
    double dz = 0.01; // z step over which DeltaM is computed
    
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz-1));
    vector<vector<vector<double> > > dMdt(Nz, vector<vector<double> > (NM, vector<double> (4,0.0)));
    
    double z = zmin;
    double Mj, Mjp1, dMdtj, dMdtjp1;
    for (int jz = 0; jz < Nz; jz++) {
        Mj = sigmalist[0][0];
        dMdtj = (1+z)*Hz(z)*DeltaM(sigmalist[0], z, dz)/dz;
        
        for (int jM = 0; jM < NM; jM++) {
            Mjp1 = sigmalist[jM+1][0];
            dMdtjp1 = (1+z)*Hz(z)*DeltaM(sigmalist[jM+1], z, dz)/dz;
            
            dMdt[jz][jM][0] = z;
            dMdt[jz][jM][1] = Mj;
            dMdt[jz][jM][2] = dMdtj;
            dMdt[jz][jM][3] = (dMdtjp1 - dMdtj)/(Mjp1 - Mj);
            
            Mj = Mjp1;
            dMdtj = dMdtjp1;
        }
        z = exp(log(z) + dlogz);
    }
    
    return dMdt;
}


// UV luminosity function, {z,M,MUV,AUV,Phi}
vector<vector<vector<double> > > cosmology::UVLFlistf(double Mt, double Mc, double epsilon, double alpha, double beta) {
    vector<vector<vector<double> > > phiUV(Nz, vector<vector<double> > (NM, vector<double> (5, 0.0)));
    double zj, Mj, dotMj, DdotMj, dndlnMj, MUVj;
    for (int jz = 0; jz < Nz; jz++) {
        zj = dotMlist[jz][0][0];
        for (int jM = 0; jM < NM; jM++) {
            Mj = dotMlist[jz][jM][1];
            dotMj = dotMlist[jz][jM][2];
            DdotMj = dotMlist[jz][jM][3];
            dndlnMj = hmflist[jz][jM][2];
            
            MUVj = MUV(Mj, dotMj, Mc, epsilon, alpha, beta);
            
            phiUV[jz][jM][0] = zj;
            phiUV[jz][jM][1] = Mj;
            phiUV[jz][jM][2] = MUVj;
            phiUV[jz][jM][3] = AUV(MUVj, zj);
            phiUV[jz][jM][4] = UVLF(Mj, dotMj, DdotMj, dndlnMj, Mt, Mc, epsilon, alpha, beta);
        }
    }
    return phiUV;
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
