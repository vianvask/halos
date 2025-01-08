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
double cosmology::sigma(double RM, double deltaH, int Nk) {
    
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


// variance of matter fluctuations, {M,sigma(M),sigma'(M)}
vector<vector<double> > cosmology::sigmalist(int Nk, int NM, double Mmin, double Mmax) {
    
    // fix deltaH to match the input sigma8
    double deltaH = sigma8/sigma(8000.0/h, 1.0, Nk);
    
    double dlogM = (log(Mmax)-log(Mmin))/(1.0*(NM-1));
    vector<vector<double> > sigma3(NM, vector<double> (3,0.0));
    double RM, M = Mmin;
    double sigman = 0.0, sigmanp = sigman;
    for (int jM = 0; jM < NM; jM++) {
        RM = pow(3.0*M/(4.0*PI*rhoM0),1.0/3.0);
        
        sigmanp = sigman;
        sigman = sigma(RM, deltaH, Nk);
        
        sigma3[jM][0] = M;
        sigma3[jM][1] = sigman;
        sigma3[jM][2] = (sigman-sigmanp)/(exp(log(M)+dlogM)-M);

        M = exp(log(M)+dlogM);
    }
    sigma3[0][2] = sigma3[1][2];
    
    return sigma3;
}


// Seth-Tormen HMF, {z,M,dn/dlnM}
vector<vector<vector<double> > > cosmology::hmflist(vector<vector<double> > &sigma3, int Nz, double zmin, double zmax) {
    function<double(double)> nuf = [&](double nu) {
        double p = 0.3;
        double q = 0.75;
        double A = 1.0/(1+(pow(2.0,-p)*tgammaf(0.5-p)/sqrt(PI)));
        return A*(1+pow(q*nu,-p))*sqrt(q*nu/(2.0*PI))*exp(-q*nu/2.0);
    };
    
    int NM = sigma3.size();
    
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz-1));
    vector<vector<vector<double> > > dndlnM(NM, vector<vector<double> > (Nz, vector<double> (3,0.0)));
    double z = zmin, M;
    for (int jz = 0; jz < Nz; jz++) {
        for (int jM = 0; jM < NM; jM++) {
            M = sigma3[jM][0];
            
            dndlnM[jz][jM][0] = z;
            dndlnM[jz][jM][1] = M;
            dndlnM[jz][jM][2] = -rhoM0*nuf(pow(deltac(z)/sigma3[jM][1],2.0))*2.0*sigma3[jM][2]/sigma3[jM][1];
        }
        z = exp(log(z) + dlogz);
    }
    
    return dndlnM;
}

// comoving distance, {z,d_c}
vector<vector<double> > cosmology::dclist(int Nz, double zmin, double zmax) {
    
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz-2));
    double z1 = 0.0, z2 = zmin;
    vector<double> d0 = {0.0, 0.0};
    
    vector<vector<double> > dlist(Nz, vector<double> (2));
    dlist[0] = d0;
    
    for (int jz = 1; jz < Nz; jz++) {
        d0[0] = z2;
        d0[1] += (z2-z1)*306.535*exp((log(1.0/Hz(z2))+log(1.0/Hz(z1)))/2.0);
        dlist[jz] = d0;
        
        z1 = z2;
        z2 = exp(log(z2)+dlogz);
    }
    return dlist;
}

// age of the Universe, {z,t}
vector<vector<double> > cosmology::tlist(int Nz, double zmin, double zmax) {
    
    double dlogz = (log(zmax)-log(zmin))/(1.0*(Nz-2));
    double z1 = 0.0, z2 = zmin;
    vector<double> d0 = {0.0, 0.0};
    
    vector<vector<double> > tlist(Nz, vector<double> (2));
    tlist[0] = d0;
    
    for (int jz = 1; jz < Nz; jz++) {
        d0[0] = z2;
        d0[1] += (z2-z1)*exp((log(1.0/((1+z2)*Hz(z2)))+log(1.0/((1+z1)*Hz(z1))))/2.0);
        
        tlist[jz] = d0;
        
        z1 = z2;
        z2 = exp(log(z2)+dlogz);
    }
    
    double tmax = tlist[Nz-1][1];
    for (int jz = 0; jz < Nz; jz++) {
        tlist[jz][1] = tmax - tlist[jz][1];
    }
    
    return tlist;
}

