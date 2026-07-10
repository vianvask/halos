#include "cosmology.h"
#include "subhalo.h"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>

array<double,2> FgNFW(double x);
double kappa0NFW(double rs, double rhos, double Sigmac);
double rmaxfNFW(cosmology &C, double zs, double zl, double M, double kappathr);
double Sigmacf(cosmology &C, double zs, double zl);

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

double safeNFWGammaCore(double x, array<double,2> &Fg) {
    if (x < 1.0e-4) return 0.5;
    return 2.0*Fg[1]/(x*x) - Fg[0];
}

// subhalo tables; evolved SHMF and normalization from Jiang & van den Bosch 2014
void Subhalo::precompute(cosmology &C, double zs, double kappathr, double kappathr_host) {
    log_Mmin = log(C.Mlist[0]);
    double log_Mmax = log(C.Mlist[C.NM-1]);
    inv_dlogM = (C.NM-1)/(log_Mmax - log_Mmin);
    
    gnorm.assign(C.Nz, vector<double>(C.NM, 0.0));
    r200h.assign(C.Nz, vector<double>(C.NM, 0.0));
    chost.assign(C.Nz, vector<double>(C.NM, 0.0));
    r_thr.assign(C.Nz, vector<double>(C.NM, 0.0));
    invRad.assign(C.Nz, vector<vector<double> >(C.NM));
    muW.assign(C.Nz, vector<vector<double> >(C.NM));
    sW.assign(C.Nz, vector<vector<double> >(C.NM));
    lyW.assign(C.Nz, vector<array<double,2> >(C.NM, {0.0, 1.0}));
    fsb.assign(C.Nz, vector<double>(C.NM, 0.0));
    
    for (int jz = 0; jz < C.Nz; jz++) {
        double zl = C.zlist[jz];
        if (zl < zs) {
            for (int jM = 0; jM < C.NM; jM++) {
                r_thr[jz][jM] = rmaxfNFW(C, zs, zl, C.Mlist[jM], kappathr);
            }
        }
    }
    
    double s = (1.0 + alpha)/omega;
    // Giocoli et al. 2007: a_f = 0.815 exp(-2 f^3)/f^0.707 with f = 1/2
    double af = 0.815*exp(-0.25)/pow(0.5, 0.707);
    double wf = sqrt(2.0*log(af + 1.0));
    double gden = gsl_sf_gamma_inc(s, beta*pow(psi_res, omega)) - gsl_sf_gamma_inc(s, beta);
    
    for (int jz = 0; jz < C.Nz; jz++) {
        double z = C.zlist[jz];
        double dcz = C.deltac(z);
        for (int jM = 0; jM < C.NM; jM++) {
            double M = C.Mlist[jM];
            if (M <= C.Mmin) continue; // no clump above the grid floor fits
            
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
            
            vector<double> NFWp = interpolate2(z, M, C.zlist, C.Mlist, C.NFWlist);
            double c = NFWp[2];
            chost[jz][jM] = c;
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

            buildWsubBin(C, zs, jz, jM, kappathr_host);
        }
    }
}

// mean and variance of unresolved clumps at fixed host impact parameter
void Subhalo::buildWsubBin(cosmology &C, double zs, int jz, int jM, double kappathr_host) {
    double zl = C.zlist[jz], M = C.Mlist[jM];
    double g = gnorm[jz][jM];
    if (g <= 0.0) return;
    double rmaxH = rmaxfNFW(C, zs, zl, M, kappathr_host);
    if (rmaxH <= 0.0) return;
    double Sigmac = Sigmacf(C, zs, zl);
    double r200 = r200h[jz][jM];
    double c = chost[jz][jM];

    double psi_min = m_floor/M;
    if (psi_min >= psi_max) return;

    double sm = (1.0 + alpha)/omega;
    double fs = g/(omega*pow(beta, sm))*(
        gsl_sf_gamma_inc(sm, beta*pow(psi_min, omega))
        - gsl_sf_gamma_inc(sm, beta*pow(psi_max, omega)));
    fsb[jz][jM] = max(0.0, min(0.95, fs));

    int NR = 128, NU = 96;
    vector<double> p2(NR, 0.0);
    double p2norm = 0.0;
    for (int i = 0; i < NR; i++) {
        double Rh = (i + 0.5)/NR;
        double umax = sqrt(max(1.0 - Rh*Rh, 0.0));
        double du = umax/(NU - 1), acc = 0.0;
        for (int k = 0; k < NU; k++) {
            double u = k*du;
            double x = sqrt(Rh*Rh + u*u);
            double B = 1.0/sqrt(pow(x/0.54, -2.5) + 1.0);
            double p3 = x*x/pow(1.0 + c*x, 2.0)*B;
            double w = (k == 0 || k == NU - 1) ? 0.5 : 1.0;
            acc += w*p3*Rh/(Rh*Rh + u*u);
        }
        p2[i] = acc*du;
        p2norm += p2[i]/NR;
    }
    if (p2norm <= 0.0) return;

    int Ny = NyW, Nd = 48, Nth = 32, Nm = 32;
    double y0 = 0.05, y1 = max(rmaxH, 2.0*y0);
    double ly0 = log(y0), dly = (log(y1) - ly0)/(Ny - 1);
    lyW[jz][jM] = {ly0, dly};
    double d0 = 0.05, d1 = rmaxH + r200;
    double dld = (log(d1) - log(d0))/(Nd - 1);
    vector<double> dg(Nd);
    for (int id = 0; id < Nd; id++) dg[id] = exp(log(d0) + id*dld);

    vector<double> fd(Ny*Nd);
    double dth = PI/(Nth - 1);
    for (int iy = 0; iy < Ny; iy++) {
        double y = exp(ly0 + iy*dly);
        for (int id = 0; id < Nd; id++) {
            double d = dg[id], ang = 0.0;
            for (int it = 0; it < Nth; it++) {
                double th = it*dth;
                double sh = sqrt(max(y*y + d*d - 2.0*y*d*cos(th), 0.0));
                double Rh = sh/r200, sig2d = 0.0;
                if (Rh < 1.0) {
                    double t = Rh*NR - 0.5, p;
                    int i = int(floor(t));
                    if (i < 0) p = p2[0]*(Rh*NR/0.5);
                    else if (i >= NR - 1) p = p2[NR - 1];
                    else p = p2[i] + (t - i)*(p2[i+1] - p2[i]);
                    p /= p2norm;
                    sig2d = p/(2.0*PI*max(Rh, 1.0e-12)*r200*r200);
                }
                double w = (it == 0 || it == Nth - 1) ? 0.5 : 1.0;
                ang += w*sig2d;
            }
            double wd = (id == 0 || id == Nd - 1) ? 0.5 : 1.0;
            fd[iy*Nd + id] = 2.0*ang*dth*d*d*dld*wd;
        }
    }

    double lp0 = log(psi_min), dlp = (log(psi_max) - lp0)/(Nm - 1);
    vector<double> J1(Nm*Ny), J2(Nm*Ny), wm(Nm), lpg(Nm), k1(Nd);
    for (int im = 0; im < Nm; im++) {
        double lp = lp0 + im*dlp;
        double psi = exp(lp), m = psi*M;
        lpg[im] = lp;
        wm[im] = g*pow(psi, alpha)*exp(-beta*pow(psi, omega));
        double rs, rhos;
        interpolateNFWMass(C, jz, m, log(m), log_Mmin, inv_dlogM, rs, rhos);
        double kappa0 = kappa0NFW(rs, rhos, Sigmac);
        for (int id = 0; id < Nd; id++)
            k1[id] = 2.0*kappa0*FgNFW(max(dg[id]/rs, 1.0e-12))[0];
        for (int iy = 0; iy < Ny; iy++) {
            double a1 = 0.0, a2 = 0.0;
            for (int id = 0; id < Nd; id++) {
                a1 += fd[iy*Nd + id]*k1[id];
                a2 += fd[iy*Nd + id]*k1[id]*k1[id];
            }
            J1[im*Ny + iy] = a1;
            J2[im*Ny + iy] = a2;
        }
    }

    muW[jz][jM].assign(Ny, 0.0);
    sW[jz][jM].assign(Ny, 0.0);
    for (int iy = 0; iy < Ny; iy++) {
        double y = exp(ly0 + iy*dly);
        int jlo = lower_bound(r_thr[jz].begin(), r_thr[jz].end(), y) - r_thr[jz].begin();
        double psi_lo = (jlo >= C.NM) ? psi_max : max(max(C.Mlist[jlo], C.Mmin), m_floor)/M;
        psi_lo = min(psi_lo, psi_max);
        if (psi_lo <= psi_min) continue;
        double lc = log(psi_lo), mu = 0.0, s2 = 0.0;
        for (int im = 0; im < Nm - 1; im++) {
            double la = lpg[im], lb = lpg[im+1];
            double F1a = wm[im]*J1[im*Ny + iy];
            double F1b = wm[im+1]*J1[(im+1)*Ny + iy];
            double F2a = wm[im]*J2[im*Ny + iy];
            double F2b = wm[im+1]*J2[(im+1)*Ny + iy];
            if (lc >= lb) {
                mu += 0.5*(F1a + F1b)*dlp;
                s2 += 0.5*(F2a + F2b)*dlp;
            } else if (lc > la) {
                double t = (lc - la)/dlp;
                mu += 0.5*(F1a + F1a + t*(F1b - F1a))*(lc - la);
                s2 += 0.5*(F2a + F2a + t*(F2b - F2a))*(lc - la);
                break;
            } else {
                break;
            }
        }
        muW[jz][jM][iy] = mu;
        sW[jz][jM][iy] = sqrt(max(s2, 0.0));
    }
}

void Subhalo::wsubTerm(int jz, int jM, double r, double &mu, double &sigma) {
    mu = 0.0;
    sigma = 0.0;
    vector<double> &m = muW[jz][jM];
    if (m.empty()) return;
    vector<double> &s = sW[jz][jM];
    double t = (log(max(r, 1.0e-12)) - lyW[jz][jM][0])/lyW[jz][jM][1];
    if (t <= 0.0) {
        mu = m.front(); sigma = s.front(); return;
    }
    if (t >= m.size() - 1) {
        mu = m.back(); sigma = s.back(); return;
    }
    int i = int(t);
    double f = t - i;
    mu = m[i] + f*(m[i+1] - m[i]);
    sigma = s[i] + f*(s[i+1] - s[i]);
}

int Subhalo::addClumps(cosmology &C, int jz, int jM, double zl, double M, double Sigmac, double r, double phi, rgen &mt, double &kappa, double &gamma1, double &gamma2) {
    double g = gnorm[jz][jM];
    if (g <= 0.0) return 0;
    
    int jlo = lower_bound(r_thr[jz].begin(), r_thr[jz].end(), r) - r_thr[jz].begin();
    if (jlo >= C.NM) return 0;
    double psi_lo = max(max(C.Mlist[jlo], C.Mmin), m_floor)/M;
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
    double invalpha = 1.0/alpha;
    
    for (int k = 0; k < Nc; k++) {
        double u = randomreal(0.0, 1.0, mt);
        double psi = pow(pa_lo + u*(pa_hi - pa_lo), invalpha);
        if (randomreal(0.0, 1.0, mt) > exp(-beta*pow(psi, omega))) continue;
        double m = psi*M;
        double logm = log(m);
        
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
        double xcl = max(d/rs, 1.0e-12);
        array<double,2> Fg = FgNFW(xcl);
        double kappac = 2.0*kappa0*Fg[0];
        double gammac = 2.0*kappa0*safeNFWGammaCore(xcl, Fg);
        
        kappa += kappac;
        gamma1 += cosphid*gammac;
        gamma2 += sinphid*gammac;
    }
    return Nc;
}
