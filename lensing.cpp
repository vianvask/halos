#include "lensing.h"

// NFW lensing amplification {lnmu, dlnmu/dlnM}
vector<double> lnmu(cosmology C, double zs, double zl, double r, double M) {
    vector<double> lnmu(2,0.0);
    
    // angular diameter distances
    double DsA = C.DL(zs)/pow(1+zs,2.0);
    double DlA = C.DL(zl)/pow(1+zl,2.0);
    double DlsA = DsA - DlA*(1+zl)/(1+zs);
    
    double Sigmac = 2.08871e16*DsA/(4.0*PI*DlA*DlsA);
    
    // NFW scale radius and density
    double rs = C.rsf(M, zl);
    double rhos = C.rhosf(M, zl);
    double Drs = C.Drsf(M, zl);
    double Drhos = C.Drhosf(M, zl);
    
    double X = r/rs;
    complex<double> X2m1 (X*X-1,0.0);
    
    double Sigma = real(2.0*rs*rhos*(sqrt(X2m1) + (PI-2.0*atan(sqrt(X2m1))) - 2.0*atan((X+1)/sqrt(sqrt(X2m1))))/pow(X2m1,3.0/2.0));
        
    double dSigmapdlnM = real(2.0*M*(Drhos*rs*X2m1 + 3.0*Drs*rhos*X*X)/pow(X2m1,2.0) + 4.0*M*(Drhos*rs*X2m1 + Drs*rhos*(4.0*X*X-1))*(PI/2.0 - atan(X2m1) - atan((1+X)/sqrt(X2m1)))/pow(X2m1,5.0/2.0));
    
    lnmu[0] = Sigma/Sigmac;
    lnmu[1] = dSigmapdlnM/Sigmac;
    
    return lnmu;
}

// TODO: probability distribution {lnmu, P1(lnmu)}
vector<vector<double> > P1(int N) {
    vector<vector<double> > P1(N, vector<double> (2, 0.0));
    
    return P1;
}
