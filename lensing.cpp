#include "lensing.h"

// NFW lensing amplification and its derivative, {lnmu, dlnmu/dlnM}
vector<double> lnmu(cosmology C, double zs, double zl, double r, double M) {
    
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
    complex<double> X2m1 (X*X-1,0.0);
    
    double Sigma = real(2.0*rs*rhos*(sqrt(X2m1) + (PI-2.0*atan(sqrt(X2m1))) - 2.0*atan((X+1)/sqrt(sqrt(X2m1))))/pow(X2m1,3.0/2.0));
    double dSigmapdlnM = real(2.0*M*(Drhos*rs*X2m1 + 3.0*Drs*rhos*X*X)/pow(X2m1,2.0) + 4.0*M*(Drhos*rs*X2m1 + Drs*rhos*(4.0*X*X-1))*(PI/2.0 - atan(X2m1) - atan((1+X)/sqrt(X2m1)))/pow(X2m1,5.0/2.0));

    vector<double> lnmu2 {Sigma/Sigmac, dSigmapdlnM/Sigmac};
    
    return lnmu2;
}


// maximal r so that lnmu(r) > lnmuthr
double rmaxf(cosmology C, double zs, double zl, double M, double lnmuthr) {
    double logrmax;
    double logr1 = log(1.0e-3), logr2 = log(1.0e5);
    
    if (lnmu(C, zs, zl, exp(logr1), M)[0] > lnmuthr) {
        while (logr2-logr1 > 0.1) {
            logrmax = (logr2+logr1)/2.0;
            if (lnmu(C, zs, zl, exp(logrmax), M)[0] > lnmuthr) {
                logr1 = logrmax;
            } else {
                logr2 = logrmax;
            }
        }
    } else {
        logrmax = logr1;
    }
    
    return exp(logrmax);
}


// number of halos withing radius rmax
double Nhf(cosmology C, double zs, double lnmuthr) {
    double Nh = 0.0;
    double rmax, M, zl, dndlnM;
    
    for (int jz = 0; jz < C.Nz; jz++) {
        for (int jM = 0; jM < C.NM; jM++) {
            zl = C.hmflist[jz][jM][0];
            if(zl < zs) {
                M = C.hmflist[jz][jM][1];
                dndlnM = C.hmflist[jz][jM][2];
                
                rmax = rmaxf(C, zs, zl, M, lnmuthr);
                
                cout << zl << "   " << M << "   " << rmax << endl;
                
                Nh += 2.0*PI*pow((1+zl)*rmax,2.0)/C.Hz(zl)*dndlnM;
            }
        }
    }
    return Nh;
}


// TODO: probability distribution {lnmu, P1(lnmu)}
vector<vector<double> > P1(cosmology C, int N, double zs, double lnmuthr) {
    vector<vector<double> > P1(N, vector<double> (2, 0.0));
    
    double Nh = Nhf(C, zs, lnmuthr);
    cout << Nh << endl;
    
    return P1;
}
