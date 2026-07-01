#include "cosmology.h"
#include "subhalo.h"
#include <gsl/gsl_sf_gamma.h>

array<double,2> FgNFW(double x);
double kappa0NFW(double rs, double rhos, double Sigmac);
double rmaxfNFW(cosmology &C, double zs, double zl, double M, double kappathr);

double linfast(double y1, double y2, double x1, double x2, double x) {
    return y1 + (x - x1)/(x2 - x1)*(y2 - y1);
}

void interpolateNFWMass(cosmology &C, int jz, double m, double logm, double logMmin, double invdlogM, double &rs, double &rhos) {
    if (m <= C.Mlist[0]) {
        rs = C.NFWlist[jz][0][0];
        rhos = C.NFWlist[jz][0][1];
    } else if (m >= C.Mlist[C.NM-1]) {
        rs = C.NFWlist[jz][C.NM-1][0];
        rhos = C.NFWlist[jz][C.NM-1][1];
    } else {
        double val = (logm - logMmin)*invdlogM;
        int jm = int(val) + 1;
        double m1 = C.Mlist[jm-1];
        double m2 = C.Mlist[jm];
        rs = linfast(C.NFWlist[jz][jm-1][0], C.NFWlist[jz][jm][0], m1, m2, m);
        rhos = linfast(C.NFWlist[jz][jm-1][1], C.NFWlist[jz][jm][1], m1, m2, m);
    }
}

// subhalo tables; evolved SHMF and normalization from Jiang & van den Bosch 2014
void Subhalo::precompute(cosmology &C, double zs, double kappathr) {
    log_Mmin = log(C.Mlist[0]);
    double log_Mmax = log(C.Mlist[C.NM-1]);
    inv_dlogM = (C.NM-1)/(log_Mmax - log_Mmin);
    
    gnorm.assign(C.Nz, vector<double>(C.NM, 0.0));
    Nsub.assign(C.Nz, vector<double>(C.NM, 0.0));
    r200h.assign(C.Nz, vector<double>(C.NM, 0.0));
    r_thr.assign(C.Nz, vector<double>(C.NM, 0.0));
    invRad.assign(C.Nz, vector<vector<double> >(C.NM));
    
    for (int jz = 0; jz < C.Nz; jz++) {
        double zl = C.zlist[jz];
        if (zl < zs) {
            for (int jM = 0; jM < C.NM; jM++) {
                r_thr[jz][jM] = rmaxfNFW(C, zs, zl, C.Mlist[jM], kappathr);
            }
        }
    }
    
    double s = (1.0 + alpha)/omega;
    double af = 0.815*exp(-1.0)/pow(0.5, 0.707);
    double wf = sqrt(2.0*log(af + 1.0));
    double gden = gsl_sf_gamma_inc(s, beta*pow(psi_res, omega)) - gsl_sf_gamma_inc(s, beta);
    
    for (int jz = 0; jz < C.Nz; jz++) {
        double z = C.zlist[jz];
        double dcz = C.deltac(z);
        for (int jM = 0; jM < C.NM; jM++) {
            double M = C.Mlist[jM];
            if (M <= 10.0*C.Mmin) continue;
            
            double sigM = interpolate(M, C.sigmalist);
            double sigH = interpolate(0.5*M, C.sigmalist);
            double dsig2 = sigH*sigH - sigM*sigM;
            if (dsig2 <= 0.0) continue;
            
            double rhs = dcz + wf*sqrt(dsig2);
            if (C.deltac(30.0) < rhs) continue;
            double zlo = z, zhi = 30.0;
            for (int it = 0; it < 60; it++) {
                double zm = 0.5*(zlo + zhi);
                if (C.deltac(zm) < rhs) {
                    zlo = zm;
                } else {
                    zhi = zm;
                }
            }
            double zf = 0.5*(zlo + zhi);
            
            double Ntau = 0.0, dz = (zf - z)/200.0;
            for (int i = 0; i < 200; i++) {
                double zz = z + (i + 0.5)*dz;
                double d = C.OmegaMz(zz) - 1.0;
                double Dvir = 18.0*PI*PI + 82.0*d - 39.0*d*d;
                Ntau += 6.006*sqrt(Dvir/178.0)/(1.0 + zz)*dz;
            }
            if (Ntau <= 0.0) continue;
            
            double fs = 0.3563/pow(Ntau, 0.6) - 0.075;
            if (fs <= 0.0 || fs >= 0.95) continue;
            double gam = omega*pow(beta, s)/gden*fs;
            gnorm[jz][jM] = gam;
            
            double psi_min = C.Mmin/M;
            if (psi_min >= psi_max) continue;
            double Nm = (gam/alpha)*(pow(psi_max, alpha) - pow(psi_min, alpha));
            if (Nm <= 0.0) continue;
            Nsub[jz][jM] = Nm;
            
            vector<double> NFWp = interpolate2(z, M, C.zlist, C.Mlist, C.NFWlist);
            double c = NFWp[2];
            r200h[jz][jM] = NFWp[0]*c;
            
            // anti-biased subhalo radial profile as in Han et al. 2016
            int Nx = 4000;
            vector<double> xs(Nx), cdf(Nx);
            double acc = 0.0;
            for (int i = 0; i < Nx; i++) {
                xs[i] = double(i)/(Nx - 1);
                double x = xs[i];
                double B = (x > 0.0) ? 1.0/sqrt(pow(x/0.54, -2.5) + 1.0) : 0.0;
                double w = x*x/pow(1.0 + c*x, 2.0)*B;
                if (i > 0) {
                    double xp = xs[i-1];
                    double Bp = (xp > 0.0) ? 1.0/sqrt(pow(xp/0.54, -2.5) + 1.0) : 0.0;
                    double wp = xp*xp/pow(1.0 + c*xp, 2.0)*Bp;
                    acc += 0.5*(w + wp)*(x - xp);
                }
                cdf[i] = acc;
            }
            
            vector<double> xu(Nu);
            for (int k = 0; k < Nu; k++) {
                double u = acc*double(k)/(Nu - 1);
                int lo = 0, hi = Nx - 1;
                while (hi - lo > 1) {
                    int mid = (lo + hi)/2;
                    if (cdf[mid] < u) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                double t = (u - cdf[lo])/(cdf[hi] - cdf[lo] + 1.0e-300);
                xu[k] = xs[lo] + t*(xs[hi] - xs[lo]);
            }
            invRad[jz][jM] = xu;
        }
    }
}

double Subhalo::resolvedFraction(cosmology &C, int jz, int jM, double M, double r) {
    double g = gnorm[jz][jM];
    if (g <= 0.0) return 0.0;
    
    int jlo = lower_bound(r_thr[jz].begin(), r_thr[jz].end(), r) - r_thr[jz].begin();
    if (jlo >= C.NM) return 0.0;
    double psi_lo = max(C.Mlist[jlo], C.Mmin)/M;
    if (psi_lo >= psi_max) return 0.0;
    
    double fs = g*(pow(psi_max, 1.0 + alpha) - pow(psi_lo, 1.0 + alpha))/(1.0 + alpha);
    return max(0.0, min(0.95, fs));
}

int Subhalo::addClumps(cosmology &C, int jz, int jM, double zl, double M, double Sigmac, double r, double phi, rgen &mt, double &kappa, double &gamma1, double &gamma2) {
    double g = gnorm[jz][jM];
    if (g <= 0.0) return 0;
    
    int jlo = lower_bound(r_thr[jz].begin(), r_thr[jz].end(), r) - r_thr[jz].begin();
    if (jlo >= C.NM) return 0;
    double psi_lo = max(C.Mlist[jlo], C.Mmin)/M;
    if (psi_lo >= psi_max) return 0;
    
    double pa_lo = pow(psi_lo, alpha);
    double pa_hi = pow(psi_max, alpha);
    double Nres = (g/alpha)*(pa_hi - pa_lo);
    if (Nres <= 0.0) return 0;
    
    poisson_distribution<int> PNsub(Nres);
    int Nc = PNsub(mt);
    if (Nc <= 0) return 0;
    
    vector<double> xcdf = invRad[jz][jM];
    double r200 = r200h[jz][jM];
    double rcos = r*cos(phi), rsin = r*sin(phi);
    double logM = log(M);
    double invalpha = 1.0/alpha;
    
    for (int k = 0; k < Nc; k++) {
        double u = randomreal(0.0, 1.0, mt);
        double logm = invalpha*log(pa_lo + u*(pa_hi - pa_lo)) + logM;
        double m = exp(logm);
        
        double ur = randomreal(0.0, 1.0, mt);
        double tt = ur*(Nu - 1);
        int i = int(tt);
        if (i >= Nu - 1) i = Nu - 2;
        double x = xcdf[i] + (tt - i)*(xcdf[i+1] - xcdf[i]);
        double r3d = x*r200;
        double cth = randomreal(-1.0, 1.0, mt);
        double psaz = randomreal(0.0, 2.0*PI, mt);
        double R2d = r3d*sqrt(1.0 - cth*cth);
        
        double dx = rcos - R2d*cos(psaz);
        double dy = rsin - R2d*sin(psaz);
        double d = sqrt(dx*dx + dy*dy);
        double invd = (d > 1.0e-30) ? 1.0/d : 0.0;
        double cosphid = dx*invd;
        double sinphid = dy*invd;
        
        double rs, rhos;
        interpolateNFWMass(C, jz, m, logm, log_Mmin, inv_dlogM, rs, rhos);
        double kappa0 = kappa0NFW(rs, rhos, Sigmac);
        double xcl = d/rs;
        array<double,2> Fg = FgNFW(xcl);
        double kappac = 2.0*kappa0*Fg[0];
        double gammac = 2.0*kappa0*(2.0*Fg[1]/(xcl*xcl) - Fg[0]);
        
        kappa += kappac;
        gamma1 += cosphid*gammac;
        gamma2 += sinphid*gammac;
    }
    return Nc;
}
