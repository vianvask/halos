#include "lensing.h"

double acot(double x) {
    return PI/2.0 - atan(x);
}

// NFW lensing amplification and its derivative, {lnmu, dlnmu/dr}
vector<double> lnmu(cosmology &C, double zs, double zl, double r, double M) {
    
    // angular diameter distances
    double DsA = C.DL(zs)/pow(1+zs,2.0);
    double DlA = C.DL(zl)/pow(1+zl,2.0);
    double DlsA = DsA - DlA*(1+zl)/(1+zs);
    
    double Sigmac = 2.08871e16*DsA/(4.0*PI*DlA*DlsA);
    
    // NFW scale radius and density
    int jz = (int) round((C.Nz-1)*log(zl/C.zmin)/log(C.zmax/C.zmin));
    int jM = (int) round((C.NM-1)*log(M/C.Mmin)/log(C.Mmax/C.Mmin));
    vector<double> zMNFW = C.NFWlist[jz][jM];
    double rs = zMNFW[2];
    double Drs = zMNFW[3];
    double rhos = zMNFW[4];
    double Drhos = zMNFW[5];

    double X = r/rs;
    double Y, Sigma, dSigmapdr;
    if (X>1) {
        Y = X*X-1;
        Sigma = 2.0*rs*rhos*(sqrt(Y) + 2.0*acot(sqrt(Y)) - 2.0*atan((X+1)/sqrt(Y)))/pow(Y,3.0/2.0);
        dSigmapdr = 2.0*rs*rhos*(2.0 - 2.0*pow(X,4.0) + Y + 6.0*sqrt(Y)*(Y+1)*(-acot(sqrt(Y)) + atan((X+1)/sqrt(Y))))/(X*pow(Y,3.0));
    } else {
        Y = 1-X*X;
        Sigma = 2.0*rs*rhos*(-sqrt(Y) + log((1+sqrt(Y))/X))/pow(Y,3.0/2.0);
        dSigmapdr = 2.0*rs*rhos*(-sqrt(Y) + log((1.0+sqrt(Y))/X))/pow(Y,3.0/2.0);
    }
    vector<double> lnmu2 {Sigma/Sigmac, dSigmapdr/Sigmac};
    
    return lnmu2;
}


// maximal r so that lnmu(r) > lnmuthr
double rmaxf(cosmology &C, double zs, double zl, double M, double lnmuthr) {
    double rmax;
    double logr1 = log(1.0e-3), logr2 = log(1.0e5);
        
    if (lnmu(C, zs, zl, exp(logr1), M)[0] > lnmuthr) {
        while (logr2-logr1 > 0.1) {
            rmax = exp((logr2+logr1)/2.0);
            if (lnmu(C, zs, zl, rmax, M)[0] > lnmuthr) {
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


// number of halos withing radius rmax
double Nhf(cosmology &C, double zs, double lnmuthr) {
    double Nh = 0.0;
    double rmax, M, zl, dndlnM;
    
    for (int jz = 0; jz < C.Nz; jz++) {
        for (int jM = 0; jM < C.NM; jM++) {
            zl = C.hmflist[jz][jM][0];
            if(zl < zs) {
                M = C.hmflist[jz][jM][1];
                dndlnM = C.hmflist[jz][jM][2];
                
                rmax = rmaxf(C, zs, zl, M, lnmuthr);
                                
                Nh += 2.0*PI*pow((1+zl)*rmax,2.0)/C.Hz(zl)*dndlnM;
            }
        }
    }
    return Nh;
}


// probability distribution normalized to N, {lnmu, dN/dlnmu}
vector<vector<double> > dNdlnmu(cosmology &C, int N, double zs, double lnmuthr, double lnmumax) {
    vector<vector<double> > N1(N, vector<double> (2, 0.0));
    
    double zl, M, dndlnM, rmax, dlnmudr, Nh;
    double dloglnmu = (log(lnmumax) - log(lnmuthr))/(1.0*N);
    double x = lnmuthr;
    
    for (int jP = 0; jP < N; jP++) {
        Nh = 0.0;
        for (int jz = 0; jz < C.Nz; jz++) {
            for (int jM = 0; jM < C.NM; jM++) {
                zl = C.hmflist[jz][jM][0];
                if(zl < zs) {
                    M = C.hmflist[jz][jM][1];
                    dndlnM = C.hmflist[jz][jM][2];
                    
                    rmax = rmaxf(C, zs, zl, M, x);
                    dlnmudr = lnmu(C, zs, zl, rmax, M)[1];
                    
                    Nh += 2.0*PI*pow(1+zl,2.0)*rmax/C.Hz(zl)*dndlnM/abs(dlnmudr);
                }
            }
        }

        N1[jP][0] = x;
        N1[jP][1] = Nh;
        
        x = exp(log(x) + dloglnmu);
    }
    return N1;
}
