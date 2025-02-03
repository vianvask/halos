#include "cosmology.h"
#include "lensing.h"

double acot(double x) {
    return PI/2.0 - atan(x);
}

// NFW lensing amplification and its derivative, {kappa, dkappa/dr}
// r_S = source size, r = impact parameter
vector<double> kappa(cosmology &C, double zs, double zl, double r, double M, double rS) {
    
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
    
    double Sigma = 0.0, dSigmapdr = 0.0;
    double X, Y;
    if (rS > 0.0) {
        // average over the source projection
        int Navg = 16;
        double rSproj = rS*DlA/DsA;
        double R, theta;
        for (int j = 0; j < Navg; j++) {
            // generate randon point inside the projection of the source in the lens plane
            R = rSproj*sqrt(randomreal(0.0,1.0));
            theta = randomreal(0.0,PI);
            
            X = sqrt(pow(r,2.0) + pow(R,2.0) - 2.0*r*R*cos(theta))/rs;
            if (X > 1.0) {
                Y = X*X-1;
                Sigma += 2*rs*rhos*(sqrt(Y) + 2*acot(sqrt(Y)) - 2*atan((X+1)/sqrt(Y)))/pow(Y,3.0/2.0);
                dSigmapdr += 2*rhos*(2 - 2*pow(X,4.0) + Y + 6*sqrt(Y)*(Y+1)*(-acot(sqrt(Y)) + atan((X+1)/sqrt(Y))))/(X*pow(Y,3.0));
            } else {
                Y = 1-X*X;
                Sigma += 2*rs*rhos*(-sqrt(Y) + log((1+sqrt(Y))/X))/pow(Y,3.0/2.0);
                dSigmapdr += 2*rhos*(-2 + 2*pow(X,4.0) + Y - 3*sqrt(Y)*(Y-1)*log((1+sqrt(Y))/X))/(X*pow(Y,3.0));
            }
        }
        Sigma = Sigma/(1.0*Navg);
        dSigmapdr = dSigmapdr/(1.0*Navg);
    } else {
        X = r/rs;
        if (X > 1.0) {
            Y = X*X-1;
            Sigma = 2*rs*rhos*(sqrt(Y) + 2*acot(sqrt(Y)) - 2*atan((X+1)/sqrt(Y)))/pow(Y,3.0/2.0);
            dSigmapdr = 2*rhos*(2 - 2*pow(X,4.0) + Y + 6*sqrt(Y)*(Y+1)*(-acot(sqrt(Y)) + atan((X+1)/sqrt(Y))))/(X*pow(Y,3.0));
        } else {
            Y = 1-X*X;
            Sigma = 2*rs*rhos*(-sqrt(Y) + log((1+sqrt(Y))/X))/pow(Y,3.0/2.0);
            dSigmapdr = 2*rhos*(-2 + 2*pow(X,4.0) + Y - 3*sqrt(Y)*(Y-1)*log((1+sqrt(Y))/X))/(X*pow(Y,3.0));
        }
    }
    
    vector<double> kappa2 {Sigma/Sigmac, dSigmapdr/Sigmac};
    
    return kappa2;
}


// maximal r so that kappa(r) > kappathr
double rmaxf(cosmology &C, double zs, double zl, double M, double kappathr, double rS) {
    double rmax;
    double logr1 = log(1.0e-3), logr2 = log(1.0e5);
        
    if (kappa(C, zs, zl, exp(logr1), M, rS)[0] > kappathr) {
        while (logr2-logr1 > 0.1) {
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


// number of halos withing radius rmax
double Nhf(cosmology &C, double zs, double kappathr, double rS) {
    double Nh = 0.0;
    double rmax, M, zl, dz, dlnM, dndlnM;
    
    for (int jz = 1; jz < C.Nz; jz++) {
        for (int jM = 1; jM < C.NM; jM++) {
            zl = C.HMFlist[jz][jM][0];
            if(zl < zs) {
                dz = zl - C.HMFlist[jz-1][jM][0];
                M = C.HMFlist[jz][jM][1];
                dlnM = log(M) - log(C.HMFlist[jz][jM-1][1]);
                dndlnM = C.HMFlist[jz][jM][2];
                
                rmax = rmaxf(C, zs, zl, M, kappathr, rS);
                Nh += 306.535*PI*pow((1+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
            }
        }
    }
    return Nh;
}


// single lens probability distribution normalized to number of haloes, {kappa, dN/dkappa}
vector<vector<double> > dNdkappa(cosmology &C, int Nx, double zs, double kappamax, double rS) {
    vector<vector<double> > N1(Nx, vector<double> (3, 0.0));
    
    // find threshold kappa that gives N_h = 50
    double kappathr;
    double kappa1 = 1.0e-6, kappa2 = 1.0;
    while (log10(kappa2)-log10(kappa1) > 0.02) {
        kappathr = pow(10.0, (log10(kappa1)+log10(kappa2))/2.0);
        if (Nhf(C, zs, kappathr, rS) > 50) {
            kappa1 = kappathr;
        } else {
            kappa2 = kappathr;
        }
    }
    
    double zl, M, dz, dlnM, dndlnM, rmax, dkappadr;
    double dlnkappa = (log(kappamax) - log(kappathr))/(1.0*(Nx-1));
    double x = kappathr, xp = x;
    
    // compute the PDF and CDF of kappa 
    double Nhcum = 0.0;
    double Nh = 0.0, Nhp = 0.0;
    for (int jx = 0; jx < Nx; jx++) {
        for (int jz = 1; jz < C.Nz; jz++) {
            for (int jM = 1; jM < C.NM; jM++) {
                zl = C.HMFlist[jz][jM][0];
                if(zl < zs) {
                    dz = zl - C.HMFlist[jz-1][jM][0];
                    M = C.HMFlist[jz][jM][1];
                    dlnM = log(M) - log(C.HMFlist[jz][jM-1][1]);
                    dndlnM = C.HMFlist[jz][jM][2];
                    
                    rmax = rmaxf(C, zs, zl, M, x, rS);
                    if (rmax > 0.0) {
                        dkappadr = kappa(C, zs, zl, rmax, M, rS)[1];
                    } else {
                        dkappadr = 1.0;
                    }
                    
                    Nh += 306.535*2.0*PI*pow(1+zl,2.0)*rmax/C.Hz(zl)/abs(dkappadr)*dndlnM*dlnM*dz;
                }
            }
        }
        
        Nhcum += exp((log(Nh) + log(Nhp))/2.0)*(x - xp);

        N1[jx][0] = x;
        N1[jx][1] = Nh;
        N1[jx][2] = Nhcum;
        
        xp = x;
        x = exp(log(x) + dlnkappa);
        
        Nhp = Nh;
        Nh = 0.0;
    }

    return N1;
}


// probability distribution of lnmu, {lnmu, dP/dlnmu}
vector<vector<double> > Plnmuf(cosmology &C, int Nx, double zs, double kappamax, double rS, int Nreal, int Nbins, rgen &mt) {

    vector<vector<double> > N1list = dNdkappa(C, Nx, zs, kappamax, rS);
    double Nh = N1list[N1list.size()-1][2];
    
    vector<vector<double> > Fm1(Nx,vector<double> (2,0.0));
    for (int j = 0; j < Nx; j++) {
        Fm1[j][0] = N1list[j][2]/Nh;
        Fm1[j][1] = log(N1list[j][0]);
    }
    
    vector<double> lnmu(Nreal, 0.0);
    vector<vector<double> > Plnmu(Nbins, vector<double> (2, 0.0));
    double x, kappatot, lnmutot, lnmumin, lnmumax, lnmumean, dlnmu, dP = 1.0/(1.0*Nreal);
    
    // generate realizations of kappa
    lnmumin = 1.0;
    lnmumax = 0.0;
    lnmumean = 0.0;
    for (int j = 0; j<Nreal; j++) {
        kappatot = 0.0;
        for (int jh = 0; jh < Nh; jh++) {
            x = randomreal(0.0,1.0,mt);
            kappatot += exp(interpolate(x,Fm1));
        }
        
        lnmutot = log(pow(1.0-kappatot,-2.0)); // mu ≈ (1-kappa)^-2 (see 1106.3823)
        if (lnmutot > lnmumax) {
            lnmumax = lnmutot;
        }
        if (lnmutot < lnmumin) {
            lnmumin = lnmutot;
        }
        lnmumean += lnmutot;
        lnmu[j] = lnmutot;
    }
    lnmumean = lnmumean/(1.0*Nreal);
    
    // binning
    dlnmu = (lnmumax-lnmumin)/(1.0*(Nbins-1));
    for (int j = 0; j < Nbins; j++) {
        Plnmu[j][0] = (lnmumin-lnmumean) + j*dlnmu;
    }
    for (int j = 0; j < Nreal; j++) {
        Plnmu[(int) round((Nbins-1)*(lnmu[j]-lnmumin)/(lnmumax-lnmumin))][1] += dP/dlnmu;
    }
    
    return Plnmu;
}


// read or generate lensing amplification distribution
vector<vector<vector<double> > > getPlnmu(cosmology &C, vector<double> &Zlist, double rS, int Nkappa, int Nreal, int Nbins) {
    vector<vector<double> > Plnmu;
    vector<vector<vector<double> > > Plnmuz;
    
    double z;
    
    ifstream infile;
    infile.open("Plnmu.dat");
    if (infile) {        
        Zlist.clear();
        vector<double> tmp(2,0.0);
        int jA = 0;
        double A;
        z = 0;
        while (infile >> A) {
            if (jA == 0) {
                if (A > z) {
                    if (z > 0) {
                        Plnmuz.push_back(Plnmu);
                        Plnmu.clear();
                    }
                    z = A;
                    Zlist.push_back(z);
                }
                jA++;
            } else {
                tmp[jA-1] = A;
                jA++;
            }
            if (jA == 3) {
                Plnmu.push_back(tmp);
                jA = 0;
            }
        }
        Plnmuz.push_back(Plnmu);
        Plnmu.clear();
    }
    else {
        ofstream outfile;
        outfile.open("Plnmu.dat");
        
        int Nkappa = 40;
        
        rgen mt(time(NULL)); // random number generator
        
        for (int jz = 0; jz < Zlist.size(); jz++) {
            z = Zlist[jz];
            cout << "z = " << z << endl;
            Plnmu = Plnmuf(C, Nkappa, z, 1.0, rS, Nreal, Nbins, mt);
            
            // output dP/dlnmu
            for (int jb = 0; jb < Nbins; jb++) {
                outfile << z << "   " << Plnmu[jb][0] << "   " << Plnmu[jb][1] << endl;
            }
            Plnmuz.push_back(Plnmu);
        }
        outfile.close();
    }
    infile.close();
    
    return Plnmuz;
}
