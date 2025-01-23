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
    if (X > 1.0) {
        Y = X*X-1;
        Sigma = 2*rs*rhos*(sqrt(Y) + 2*acot(sqrt(Y)) - 2*atan((X+1)/sqrt(Y)))/pow(Y,3.0/2.0);
        dSigmapdr = 2*rhos*(2 - 2*pow(X,4.0) + Y + 6*sqrt(Y)*(Y+1)*(-acot(sqrt(Y)) + atan((X+1)/sqrt(Y))))/(X*pow(Y,3.0));
    } else {
        Y = 1-X*X;
        Sigma = 2*rs*rhos*(-sqrt(Y) + log((1+sqrt(Y))/X))/pow(Y,3.0/2.0);
        dSigmapdr = 2*rhos*(-2 + 2*pow(X,4.0) + Y - 3*sqrt(Y)*(Y-1)*log((1+sqrt(Y))/X))/(X*pow(Y,3.0));
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
    double rmax, M, zl, dz, dlnM, dndlnM;
    
    for (int jz = 1; jz < C.Nz; jz++) {
        for (int jM = 1; jM < C.NM; jM++) {
            zl = C.hmflist[jz][jM][0];
            if(zl < zs) {
                dz = zl - C.hmflist[jz-1][jM][0];
                M = C.hmflist[jz][jM][1];
                dlnM = log(M) - log(C.hmflist[jz][jM-1][1]);
                dndlnM = C.hmflist[jz][jM][2];
                
                rmax = rmaxf(C, zs, zl, M, lnmuthr);
                Nh += 306.535*PI*pow((1+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
            }
        }
    }
    return Nh;
}


// probability distribution normalized to N, {lnmu, dN/dlnmu}
vector<vector<double> > dNdlnmu(cosmology &C, int Nx, double zs, double lnmumax) {
    vector<vector<double> > N1(Nx, vector<double> (3, 0.0));
    
    // find threshold lnmu that gives N_h = 50
    double lnmuthr;
    double lnmu1 = 1.0e-6, lnmu2 = 1.0;
    while (log10(lnmu2)-log10(lnmu1) > 0.02) {
        lnmuthr = pow(10.0, (log10(lnmu1)+log10(lnmu2))/2.0);
        if (Nhf(C, zs, lnmuthr) > 50) {
            lnmu1 = lnmuthr;
        } else {
            lnmu2 = lnmuthr;
        }
    }
    
    double zl, M, dz, dlnM, dndlnM, rmax, dlnmudr;
    double dloglnmu = (log(lnmumax) - log(lnmuthr))/(1.0*(Nx-1));
    double x = lnmuthr, xp = x;
    
    // compute the PDF and CDF of lnmu 
    double Nhcum = 0.0;
    double Nh = 0.0, Nhp = 0.0;
    for (int jx = 0; jx < Nx; jx++) {
        for (int jz = 1; jz < C.Nz; jz++) {
            for (int jM = 1; jM < C.NM; jM++) {
                zl = C.hmflist[jz][jM][0];
                if(zl < zs) {
                    dz = zl - C.hmflist[jz-1][jM][0];
                    M = C.hmflist[jz][jM][1];
                    dlnM = log(M) - log(C.hmflist[jz][jM-1][1]);
                    dndlnM = C.hmflist[jz][jM][2];
                    
                    rmax = rmaxf(C, zs, zl, M, x);
                    if (rmax > 0.0) {
                        dlnmudr = lnmu(C, zs, zl, rmax, M)[1];
                    } else {
                        dlnmudr = 1.0;
                    }
                    
                    Nh += 306.535*2.0*PI*pow(1+zl,2.0)*rmax/C.Hz(zl)/abs(dlnmudr)*dndlnM*dlnM*dz;
                }
            }
        }
        
        Nhcum += exp((log(Nh) + log(Nhp))/2.0)*(x - xp);

        N1[jx][0] = x;
        N1[jx][1] = Nh;
        N1[jx][2] = Nhcum;
        
        xp = x;
        x = exp(log(x) + dloglnmu);
        
        Nhp = Nh;
        Nh = 0.0;
    }

    return N1;
}
