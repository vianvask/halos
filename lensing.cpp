#include "lensing.h"

// TODO: NFW scale radius and density
double rsf(double M, double z) {
    double r;
    
    return r;
}
double rhosf(double M, double z) {
    double rho;
    
    return rho;
}
double Drsf(double M, double z) {
    double Dr;
    
    return Dr;
}
double Drhosf(double M, double z) {
    double Drho;
    
    return Drho;
}

// NFW lensing amplification {lnmu, dlnmu/dlnM}
vector<double> lnmu(cosmology C, double zs, double zl, double r, double M) {
    vector<double> lnmu(2,0.0);
    
    double Sigmac = C.DL(zs)/(4.0*PI*C.DL(zl)*(C.DL(zs) - C.DL(zl)*(1+zl)/(1+zs)));
    
    double rs = rsf(M, zl);
    double rhos = rhosf(M, zl);
    double Drs = Drsf(M, zl);
    double Drhos = Drhosf(M, zl);
    
    double X = r/rs;
    double X2m1 = X*X-1;
    
    double Sigma = 2.0*rs*rhos*(sqrt(X2m1) + (PI-2.0*atan(sqrt(X2m1))) - 2.0*atan((X+1)/sqrt(sqrt(X2m1))))/pow(X2m1,3.0/2.0);
    
    double dSigmapdlnM = 2.0*M*(Drhos*rs*X2m1 + 3.0*Drs*rhos*X*X)/pow(X2m1,2.0) + 4.0*M*(Drhos*rs*X2m1 + Drs*rhos*(4.0*X*X-1))*(PI/2.0 - atan(X2m1) - atan((1+X)/sqrt(X2m1)))/pow(X2m1,5.0/2.0);
    
    lnmu[0] = Sigma/Sigmac;
    lnmu[1] = dSigmapdlnM/Sigmac;
    
    return lnmu;
}

// TODO: probability distribution {lnmu, P1(lnmu)}
vector<vector<double> > P1(int N) {
    vector<vector<double> > P1(N, vector<double> (2, 0.0));
    
    return P1;
}
