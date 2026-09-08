#include "cosmology.h"
#include "subhalo.h"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>

array<double,2> FgNFW(double x);
double kappa0NFW(double rs, double rhos, double Sigmac);
double rmaxfNFW(cosmology &C, double zs, double zl, double M, double kappathr);
double Sigmacf(cosmology &C, double zs, double zl);

// anti-biased radial bias B(x) = [1 + (x/x0)^-p]^-1/2, calibrated in virial
// radius units and converted to r200 units per host by biasScaleR200
static const double BIAS_X0_RVIR = 0.86; // transition scale in r_vir units
static const double BIAS_EXP = 2.5; // p

// Green et al. 2021 define B as a ratio of volume number densities, so B
// multiplies rho_NFW and the shell volume cancels one power of x against the
// NFW cusp: dN/dx = 4 pi r^2 rho_NFW B ~ x B/(1+cx)^2, not x^2 B/(1+cx)^2.
// Moves clumps inward at fixed total mass per host, so it is not bitwise.

double linfast(double y1, double y2, double x1, double x2, double x) {
    return y1 + (x - x1)/(x2 - x1)*(y2 - y1);
}

// r_vir/r200 = c_vir/c200 for an NFW halo, both radii sharing r_s. c_vir solves
// 200 c200^3/m(c200) = Dvir c_vir^3/m(c_vir) with Dvir from Bryan & Norman 1998.
// Dvir < 200 in LCDM, so eta > 1 and the bracket (c200, 3 c200] always holds.
double nfwMu(double y) { return log(1.0 + y) - y/(1.0 + y); }

double etaVirTo200(cosmology &C, double c200, double z) {
    const double d = C.OmegaMz(z) - 1.0;
    const double Dvir = 18.0*PI*PI + 82.0*d - 39.0*d*d;
    const double target = 200.0*c200*c200*c200/nfwMu(c200)/Dvir;
    double lo = c200, hi = 3.0*c200; // c^3/m(c) is monotone increasing
    for (int it = 0; it < 80; it++) {
        double mid = 0.5*(lo + hi);
        if (mid*mid*mid/nfwMu(mid) < target) lo = mid; else hi = mid;
    }
    return 0.5*(lo + hi)/c200;
}

// transition scale of B(x) in x = r/r200; the fit value is in r_vir units
double biasScaleR200(cosmology &C, double c200, double z) {
    return BIAS_X0_RVIR*etaVirTo200(C, c200, z);
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
void Subhalo::precompute(cosmology &C, double zs, double kappathr, double kappathr_host,
                         bool build_restricted) {
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
    // virial convention tables; the legacy values 1.0 and M make the
    // expressions below reduce exactly to the pre-virial form
    xmaxh.assign(C.Nz, vector<double>(C.NM, 1.0));
    Mpsih.assign(C.Nz, vector<double>(C.NM, 0.0));
    Mgrid_.assign(C.Nz, vector<double>(C.NM, 0.0));
    // restricted sampling tables, subhalo_model 5
    p2p.assign(C.Nz, vector<vector<double> >(C.NM));
    propCum.assign(C.Nz, vector<vector<double> >(C.NM));
    sig2dmax.assign(C.Nz, vector<double>(C.NM, 0.0));
    lpq.assign(C.Nz, vector<array<double,2> >(C.NM, {0.0, 1.0}));
    restricted_built = false;

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

            // virial scales: M_vir/M200 = mu(c_vir)/mu(c200), r_vir/r200 = eta
            Mgrid_[jz][jM] = M;
            if (virial) {
                const double eta = etaVirTo200(C, c, z);
                xmaxh[jz][jM] = eta;
                Mpsih[jz][jM] = M*nfwMu(eta*c)/nfwMu(c);
            } else {
                Mpsih[jz][jM] = M; // psi = m/M200, extent x <= 1
            }
            const double xmax = xmaxh[jz][jM];

            // anti-biased subhalo radial profile as in Han et al. 2016, sampled
            // over x in [0, xmax]: r200 in legacy mode, r_vir in virial mode
            const double x0 = biasScaleR200(C, c, z);
            int Nx = 4000;
            vector<double> xs(Nx), cdf(Nx);
            double acc = 0.0;
            for (int i = 0; i < Nx; i++) {
                xs[i] = double(i)/(Nx - 1)*xmax;
                double x = xs[i];
                double B = (x > 0.0) ? 1.0/sqrt(pow(x/x0, -BIAS_EXP) + 1.0) : 0.0;
                double w = x/pow(1.0 + c*x, 2.0)*B;
                if (i > 0) {
                    double xp = xs[i-1];
                    double Bp = (xp > 0.0) ? 1.0/sqrt(pow(xp/x0, -BIAS_EXP) + 1.0) : 0.0;
                    double wp = xp/pow(1.0 + c*xp, 2.0)*Bp;
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

            // only model 3 has an unresolved band and needs the Wsub tables
            if (kappathr_host > 0.0) {
                buildWsubBin(C, zs, jz, jM, kappathr_host);
            }
            if (build_restricted) {
                buildRestrictedBin(C, jz, jM);
            }
        }
    }
    restricted_built = build_restricted;
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

    double psi_min = (psi_min_fixed > 0.0) ? psi_min_fixed : (m_floor/M);
    if (psi_min >= psi_max) return;

    double sm = (1.0 + alpha)/omega;
    double fs = g/(omega*pow(beta, sm))*(
        gsl_sf_gamma_inc(sm, beta*pow(psi_min, omega))
        - gsl_sf_gamma_inc(sm, beta*pow(psi_max, omega)));
    fsb[jz][jM] = max(0.0, min(0.95, fs));

    // projected surface density on a cell-centered Rhat grid. Model 3 is legacy
    // convention only, virial being gated to models 4 and 5, so the extent is 1
    int NR = 128, NU = 96;
    const double x0 = biasScaleR200(C, c, zl);
    vector<double> p2(NR, 0.0);
    double p2norm = 0.0;
    for (int i = 0; i < NR; i++) {
        double Rh = (i + 0.5)/NR;
        double umax = sqrt(max(1.0 - Rh*Rh, 0.0));
        double du = umax/(NU - 1), acc = 0.0;
        for (int k = 0; k < NU; k++) {
            double u = k*du;
            double x = sqrt(Rh*Rh + u*u);
            double B = 1.0/sqrt(pow(x/x0, -BIAS_EXP) + 1.0);
            double p3 = x/pow(1.0 + c*x, 2.0)*B;
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
        // same dynamic floor as addClumps and unresolvedMass; the three must agree
        double psi_lo = (jlo >= C.NM) ? psi_max : max(C.Mlist[jlo], C.Mmin)/M;
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

// mean unresolved bound mass at ray distance r, for the mass conserving carve.
// The carved host M - sum_i m_i - M_u(r) reproduces the mean host mass
// (1-f_s,b)M while conserving the total halo mass M in every realization.
// Model 3 only: brute and model 4 have every clump explicit, and model 5 keeps
// the sub-threshold mass in the smooth host by construction.
double Subhalo::unresolvedMass(cosmology &C, int jz, int jM, double r, double M) const {
    double g = gnorm[jz][jM];
    if (g <= 0.0) return 0.0;
    double f_b = fsb[jz][jM]; // full bound fraction over [psi_min, psi_max]
    if (f_b <= 0.0) return 0.0;
    // smallest clump resolved at this distance; below it everything is unresolved
    const vector<double> &rth = r_thr[jz];
    int jlo = int(lower_bound(rth.begin(), rth.end(), r) - rth.begin());
    double psi_lo = (jlo >= C.NM) ? psi_max : max(C.Mlist[jlo], C.Mmin)/M;
    psi_lo = min(psi_lo, psi_max);
    const double psi_min = (psi_min_fixed > 0.0) ? psi_min_fixed : (m_floor/M);
    if (psi_lo <= psi_min) return 0.0; // everything resolved
    // resolved bound fraction above psi_lo, so M_u is the complementary mass
    const double s_m = (1.0 + alpha)/omega;
    double f_s_res = g/(omega*pow(beta, s_m))*(
        gsl_sf_gamma_inc(s_m, beta*pow(psi_lo, omega))
        - gsl_sf_gamma_inc(s_m, beta*pow(psi_max, omega)));
    f_s_res = max(0.0, min(0.95, f_s_res));
    double f_unr = f_b - f_s_res;
    return (f_unr > 0.0) ? f_unr*M : 0.0;
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                        restricted intensity clump sampling, subhalo_model 5                                                    */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// normalized projected clump number density at lens plane separation s from the
// host centre, zero outside the population extent
double Subhalo::sigma2Dclump(int jz, int jM, double s) const {
    const vector<double> &p2 = p2p[jz][jM];
    if (p2.empty()) return 0.0;
    const double r200 = r200h[jz][jM];
    const double xmax = xmaxh[jz][jM];
    const double Rh = s/r200;
    if (Rh >= xmax) return 0.0;
    const int NR = int(p2.size());
    double t = Rh/xmax*NR - 0.5; // p2 spans Rhat in [0, xmax] over NR cells
    int i = int(floor(t));
    double p;
    if (i < 0) p = p2[0]*(Rh/xmax*NR/0.5); // p2 ~ Rhat near 0
    else if (i >= NR - 1) p = p2[NR - 1];
    else p = p2[i] + (t - i)*(p2[i+1] - p2[i]);
    return p/(2.0*PI*max(Rh, 1.0e-12)*r200*r200);
}

// clump reach D(m) from r_thr, which for model 5 is built at kappathr_sub
double Subhalo::clumpReach(cosmology &C, int jz, double m) const {
    const vector<double> &rt = r_thr[jz];
    if (rt.empty()) return 0.0;
    double t = (log(m) - log_Mmin)*inv_dlogM;
    if (t <= 0.0) return rt[0];
    if (t >= C.NM - 1) return rt[C.NM - 1];
    int i = int(t);
    return rt[i] + (t - i)*(rt[i+1] - rt[i]);
}

void Subhalo::buildRestrictedBin(cosmology &C, int jz, int jM) {
    const double g = gnorm[jz][jM];
    if (g <= 0.0) return;
    const double r200 = r200h[jz][jM];
    const double c = chost[jz][jM];
    if (r200 <= 0.0) return;

    const double Mpsi = Mpsih[jz][jM]; // M200 legacy, M_vir in virial mode
    const double xmax = xmaxh[jz][jM]; // 1 legacy, eta in virial mode
    const double psi_min = (psi_min_fixed > 0.0) ? psi_min_fixed : (m_floor/Mpsi);
    if (psi_min >= psi_max) return;

    // projected clump surface density, built as in buildWsubBin but kept
    // separate so that model 5 does not depend on model 3
    const int NR = NRp, NUa = 96;
    const double x0 = biasScaleR200(C, c, C.zlist[jz]);
    vector<double> p2(NR, 0.0);
    double p2norm = 0.0;
    for (int i = 0; i < NR; i++) {
        double Rh = (i + 0.5)/NR*xmax;
        double umax = sqrt(max(xmax*xmax - Rh*Rh, 0.0));
        double du = umax/(NUa - 1), acc2 = 0.0;
        for (int k = 0; k < NUa; k++) {
            double u = k*du;
            double w = (k == 0 || k == NUa - 1) ? 0.5 : 1.0;
            double x = sqrt(Rh*Rh + u*u);
            double p3 = 0.0;
            if (x > 0.0 && x <= xmax) {
                double B = 1.0/sqrt(pow(x/x0, -BIAS_EXP) + 1.0);
                p3 = x/pow(1.0 + c*x, 2.0)*B;
            }
            acc2 += w*p3*Rh/(Rh*Rh + u*u);
        }
        p2[i] = acc2*du;
        p2norm += p2[i]*xmax/NR;
    }
    if (p2norm <= 0.0) return;
    for (int i = 0; i < NR; i++) p2[i] /= p2norm; // int_0^xmax p2 dRhat = 1
    p2p[jz][jM] = p2;

    // thinning envelope, the maximum of Sigma_n over the plane. The R -> 0 limit
    // is finite since p2 ~ Rhat, so include it alongside the cell centres
    double smax = p2[0]*NR/(0.5*xmax*2.0*PI*r200*r200);
    for (int i = 0; i < NR; i++) {
        double Rh = (i + 0.5)/NR*xmax;
        smax = max(smax, p2[i]/(2.0*PI*Rh*r200*r200));
    }
    if (!(smax > 0.0)) return;
    sig2dmax[jz][jM] = smax;

    // proposal intensity lambda(psi) = g psi^alpha min(1, pi D^2 Smax) on a
    // log-psi grid. Independent of the ray position because the envelope is
    // global, so it precomputes once per bin; the SHMF cutoff is thinned later
    const double lp0 = log(psi_min);
    const double dlp = (log(psi_max) - lp0)/(Nq - 1);
    lpq[jz][jM] = {lp0, dlp};
    vector<double> cum(Nq, 0.0);
    double acc = 0.0, prev = 0.0;
    for (int q = 0; q < Nq; q++) {
        double psi = exp(lp0 + q*dlp);
        double D = clumpReach(C, jz, psi*Mpsi);
        double pcap = min(1.0, PI*D*D*smax);
        double lam = g*pow(psi, alpha)*pcap;
        if (q > 0) acc += 0.5*(lam + prev)*dlp; // trapezoid in dlnpsi
        prev = lam;
        cum[q] = acc;
    }
    propCum[jz][jM] = cum;
}

int Subhalo::addClumpsRestricted(cosmology &C, int jz, int jM, double M, double Sigmac,
                                 double r, double phi, rgen &mt,
                                 double &kappa, double &gamma1, double &gamma2,
                                 double *mass_out) {
    const vector<double> &cum = propCum[jz][jM];
    if (cum.empty()) return 0;
    const double Ntot = cum.back();
    if (!(Ntot > 0.0)) return 0;

    poisson_distribution<int> pois(Ntot);
    int Nprop = pois(mt);
    if (Nprop <= 0) return 0;

    const double lp0 = lpq[jz][jM][0], dlp = lpq[jz][jM][1];
    const double smax = sig2dmax[jz][jM];
    const double r200 = r200h[jz][jM];
    const vector<double> &xcdf = invRad[jz][jM];
    const double rcos = r*cos(phi), rsin = r*sin(phi);

    int Nrendered = 0;
    for (int k = 0; k < Nprop; k++) {
        // clump mass from the proposal intensity by inverse CDF
        double u = randomreal(0.0, 1.0, mt)*Ntot;
        int q = int(lower_bound(cum.begin(), cum.end(), u) - cum.begin());
        if (q <= 0) q = 1;
        if (q >= Nq) q = Nq - 1;
        double f = (cum[q] > cum[q-1]) ? (u - cum[q-1])/(cum[q] - cum[q-1]) : 0.0;
        double psi = exp(lp0 + (q - 1 + f)*dlp);
        if (psi >= psi_max) continue;

        // exponential cutoff of the SHMF, by thinning as in model 4
        if (randomreal(0.0, 1.0, mt) > exp(-beta*pow(psi, omega))) continue;

        const double m = psi*Mpsih[jz][jM];
        const double D = clumpReach(C, jz, m);
        if (!(D > 0.0)) continue;

        double d, dx, dy;
        if (PI*D*D*smax < 1.0) {
            // small target: uniform in the retention disc, thinned by Sigma_n/Smax
            double ud = randomreal(0.0, 1.0, mt);
            d = D*sqrt(ud);
            double th = randomreal(0.0, 2.0*PI, mt);
            dx = d*cos(th);
            dy = d*sin(th);
            double Rx = rcos - dx, Ry = rsin - dy; // clump radius from host centre
            double Rc = sqrt(Rx*Rx + Ry*Ry);
            if (randomreal(0.0, 1.0, mt)*smax > sigma2Dclump(jz, jM, Rc)) continue;
        } else {
            // big reach: full radial profile as in model 4, then test d <= D
            double ur = randomreal(0.0, 1.0, mt);
            double tt = ur*(Nu - 1);
            int i = int(tt); if (i >= Nu - 1) i = Nu - 2;
            double x = xcdf[i] + (tt - i)*(xcdf[i+1] - xcdf[i]);
            double r3d = x*r200;
            double cth = randomreal(-1.0, 1.0, mt);
            double psaz = randomreal(0.0, 2.0*PI, mt);
            double R2d = r3d*sqrt(1.0 - cth*cth);
            dx = rcos - R2d*cos(psaz);
            dy = rsin - R2d*sin(psaz);
            d = sqrt(dx*dx + dy*dy);
            if (d > D) continue;
        }

        if (mass_out) *mass_out += m; // retained clump mass only
        double invd = (d > 1.0e-30) ? 1.0/d : 0.0;
        double cosphid = dx*invd, sinphid = dy*invd;
        double rs, rhos;
        interpolateNFWMass(C, jz, m, log(m), log_Mmin, inv_dlogM, rs, rhos);
        double kappa0 = kappa0NFW(rs, rhos, Sigmac);
        double xcl = max(d/rs, 1.0e-12);
        array<double,2> Fg = FgNFW(xcl);
        kappa += 2.0*kappa0*Fg[0];
        double gammac = 2.0*kappa0*safeNFWGammaCore(xcl, Fg);
        gamma1 += cosphid*gammac;
        gamma2 += sinphid*gammac;
        Nrendered++;
    }
    return Nrendered;
}

int Subhalo::addClumps(cosmology &C, int jz, int jM, double zl, double M, double Sigmac, double r, double phi, rgen &mt, double &kappa, double &gamma1, double &gamma2,
                       int subhalo_model, bool subhalo_brute, double *mass_out) {
    double g = gnorm[jz][jM];
    if (g <= 0.0) return 0;

    // mass scale for psi: M200 in legacy mode, M_vir in virial mode
    const double Mpsi = Mpsih[jz][jM];

    double psi_lo = 0.0;
    if (subhalo_brute || subhalo_model == 4) {
        // resolve every subhalo down to the absolute floor, no dynamic floor
        psi_lo = (psi_min_fixed > 0.0) ? psi_min_fixed : (m_floor/Mpsi);
    } else {
        // dynamic floor, the smallest clump whose reach r_thr(m) >= r, keyed to
        // the host centre distance so that it matches the host reduction
        int jlo = int(lower_bound(r_thr[jz].begin(), r_thr[jz].end(), r) - r_thr[jz].begin());
        if (jlo >= C.NM) return 0; // nothing reaches kappathr at distance r
        psi_lo = max(C.Mlist[jlo], C.Mmin)/Mpsi;
    }
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
        double m = psi*Mpsi;
        double logm = log(m);
        if (mass_out) *mass_out += m; // realized clump mass carved from the host

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
