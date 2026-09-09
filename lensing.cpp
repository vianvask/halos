#include "cosmology.h"
#include "lensing.h"
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>


double Sigmacf(cosmology &C, double zs, double zl) {
    // angular diameter distances
    double DsA = C.DL(zs)/pow(1+zs,2.0);
    double DlA = C.DL(zl)/pow(1+zl,2.0);
    double DlsA = DsA - DlA*(1+zl)/(1+zs);
    
    // Sigma_c
    return 2.08871e16*DsA/(4.0*PI*DlA*DlsA);
}


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                              NFW halos                                                                         */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// from astro-ph/0608153
array<double,2> FgNFW(double x) {
    if (x > 1) {
        double t = atan(sqrt((x-1)/(1+x))) / sqrt(x*x - 1);
        return {(1 - 2*t) / (x*x - 1), 2*t + log(x/2)};
    }
    if (x < 1) {
        double t = atanh(sqrt((1-x)/(1+x))) / sqrt(1 - x*x);
        return {(1 - 2*t) / (x*x - 1), 2*t + log(x/2)};
    }
    return {1.0/3.0, 1 + log(0.5)};
}

double kappa0NFW(double rs, double rhos, double Sigmac) {
    return rs*rhos/Sigmac;
}

// fix ellipticity using the fit of astro-ph/0508497
double epsilonNFW(cosmology &C, double z, double M) {
    double s = 0.54*pow(M/exp(interpolate(z, C.logMcharlist)), -0.05);
    return max(0.0, (1.0-s)/(1.0+s));
}

// kappa and gamma for pseudo elliptical NFW halos
array<double,2> kappagammaNFWeps(double epsilon, double kappa0, double x, double phi) {
    double a1eps = 1.0-epsilon;
    double a2eps = 1.0+epsilon;
    double x1eps = sqrt(a1eps)*cos(phi)*x;
    double x2eps = sqrt(a2eps)*sin(phi)*x;
    double xeps = sqrt(pow(x1eps,2.0) + pow(x2eps,2.0));
    double phieps = atan(x2eps/x1eps);
    
    auto Fg = FgNFW(xeps);
    double kappaeps0 = 2.0*kappa0*Fg[0];
    double gammaeps0 = 2.0*kappa0*(2.0*Fg[1]/(xeps*xeps) - Fg[0]);

    double kappaeps = kappaeps0 + epsilon*cos(2.0*phieps)*gammaeps0;
    double gammaeps = sqrt(pow(gammaeps0,2.0) + 2.0*epsilon*cos(2.0*phieps)*gammaeps0*kappaeps0 + pow(epsilon,2.0)*(pow(kappaeps0,2.0) - pow(cos(2.0*phieps)*gammaeps0,2.0)));
    
    return {kappaeps, gammaeps};
}
array<double,2> kappagammaNFW(cosmology &C, double zs, double zl, double r, double M, double phi, double epsilon) {
    double Sigmac = Sigmacf(C, zs, zl);
    
    // NFW scale radius and density
    vector<double> NFWparams = interpolate2(zl, M, C.zlist, C.Mlist, C.NFWlist);
    double rs = NFWparams[0];
    double rhos = NFWparams[1];
    
    double kappa0 = kappa0NFW(rs, rhos, Sigmac);
    
    return kappagammaNFWeps(epsilon, kappa0, r/rs, phi);
}

// maximal r so that kappa_NFW > kappa_thr
double rmaxfNFW(cosmology &C, double zs, double zl, double M, double kappathr) {
    double rmax;
    double logr1 = log(1.0e-6), logr2 = log(1.0e6);
    if (kappagammaNFW(C, zs, zl, exp(logr1), M, 0.0, 0.0)[0] > kappathr) {
        while (logr2-logr1 > 0.02) {
            rmax = exp((logr2+logr1)/2.0);
            if (kappagammaNFW(C, zs, zl, rmax, M, 0.0, 0.0)[0] > kappathr) {
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

// number of halos with kappa_NFW > kappa_thr in each z and M bin
vector<vector<vector<double> > > deltaNhfNFW(cosmology &C, double zs, double kappathr) {
    vector<vector<vector<double> > > dNh(C.Nz, vector<vector<double> > (C.NM, vector<double> (3, 0.0)));
    double zl, dz, M, Mb, dlnM, dndlnM, rmax, sigma, sigmab;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                sigma = C.sigmalist[jM][1];
                
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                rmax = rmaxfNFW(C, zs, zl, M, kappathr);
                dNh[jz][jM][0] = CLIGHT*PI*pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
                
                Mb = 2.0*PI*pow(rmax,2.0)*(C.dc(zl)-C.dc(zl-dz))*C.rhoM0;
                sigmab = interpolate(Mb, C.sigmalist);
                
                // see Baumann (5.129)
                dNh[jz][jM][1] = C.Dg(zl)*sigmab*C.halobias(zl, sigma);
                
                dNh[jz][jM][2] = rmax;
            }
        }
    }
    return dNh;
}

// total number of halos with kappa_NFW > kappa_thr
double NhfNFW(cosmology &C, double zs, double kappathr) {
    double Nh = 0.0;
    double zl, dz, M, dlnM, dndlnM, rmax;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                rmax = rmaxfNFW(C, zs, zl, M, kappathr);
                Nh += CLIGHT*PI*pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
            }
        }
    }
    return Nh;
}

// variance of kappa from weak lenses
double sigmakappaW(cosmology &C, double zs, double kappathr) {
    double kappa2 = 0.0;
    double zl, dz, M, dlnM, dndlnM, r, kappar;
    int Nr = 100;
    double dlnr = 0.01;
    double Edlnr = exp(dlnr);
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];
                                
                r = rmaxfNFW(C, zs, zl, M, kappathr);
                if (r == 0.0) {
                    r = 1.0e-6;
                }

                double SigmacW = Sigmacf(C, zs, zl);
                vector<double> NFWpW = interpolate2(zl, M, C.zlist, C.Mlist, C.NFWlist);
                double rsW = NFWpW[0];
                double kappa0W = kappa0NFW(rsW, NFWpW[1], SigmacW);

                kappar = kappathr;
                while (kappar > 0.001*kappathr) {
                    kappar = kappagammaNFWeps(0.0, kappa0W, r/rsW, 0.0)[0];
                    // Campbell's theorem for Poisson-distributed halo counts: Var = int n kappa^2,
                    // with log-annulus area element d(pi r^2) = 2 pi r^2 dlnr
                    kappa2 += CLIGHT*2.0*PI*pow((1.0+zl)*r,2.0)/C.Hz(zl)*dndlnM*pow(kappar,2.0)*dlnr*dlnM*dz;
                    r = r*Edlnr;
                }
            }
        }
    }
    return sqrt(kappa2);
}


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                           Cylindrical filaments                                                                */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// axis in the lens plane
double kappa0CYL(double rs, double rhos, double Sigmac) {
    return rs*rhos/Sigmac;
}

// uniform density inside the cylinder
double kappaCYL1(double r, double rs, double kappa0) {
    if (r > rs) {
        return 0.0;
    }
    return kappa0*2.0*sqrt(1.0-pow(r/rs,2.0));
}
double gammaCYL1(double r, double rs, double kappa0) {
    return kappaCYL1(r, rs, kappa0);
}

// 1/(1+(r/r_s)^2) profile inside the cylinder
double kappaCYL2(double r, double rs, double kappa0) {
    if (r > rs) {
        return 0.0;
    }
    return kappa0*2.0*PI/sqrt(1.0+pow(r/rs,2.0));
}
double gammaCYL2(double r, double rs, double kappa0) {
    return kappa0*4.0*PI*(sqrt(1.0+pow(r/rs,2.0)) - 1.0)*pow(rs/r,2.0) - kappaCYL2(r, rs, kappa0);
}

array<double,2> kappagammaCYL(cosmology &C, double zs, double zl, double r, double M, double phi) {
    double Sigmac = Sigmacf(C, zs, zl);
    
    // cylinder radius and length, M = mass inside radius r_s (see astro-ph/0406665)
    double rs = 1000.0*pow(M/1.0e14, 1.0/3.0);
    double L = 20000.0*pow(M/1.0e14, 1.0/3.0);
    
    // average density inside radius r_s = 10*rho_c, rough approximation of rho_s density when it is not in the lens plane
    double rhos = 14.4*C.rhoc*L/max(2.0*rs, abs(cos(phi))*L);
    
    double kappa0 = kappa0CYL(rs, rhos, Sigmac);
    
    return {kappaCYL2(r, rs, kappa0), gammaCYL2(r, rs, kappa0)};
}

// maximal r so that kappa_CYL > kappa_thr
double rmaxfCYL(cosmology &C, double zs, double zl, double M, double kappathr) {
    double rmax;
    double logr1 = log(1.0e-6), logr2 = log(1.0e6);
    if (kappagammaCYL(C, zs, zl, exp(logr1), M, 0.0)[0] > kappathr) {
        while (logr2-logr1 > 0.02) {
            rmax = exp((logr2+logr1)/2.0);
            if (kappagammaCYL(C, zs, zl, rmax, M, 0.0)[0] > kappathr) {
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

// number of halos with kappa_NFW > kappa_thr in each z and M bin
vector<vector<vector<double> > > deltaNhfCYL(cosmology &C, double zs, double kappathr) {
    vector<vector<vector<double> > > dNh(C.Nz, vector<vector<double> > (C.NM, vector<double> (3, 0.0)));
    double zl, dz, M, Mb, dlnM, dndlnM, rmax, sigma, sigmab;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                sigma = C.sigmalist[jM][1];
                
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.FMFlist[jz][jM];
                rmax = rmaxfCYL(C, zs, zl, M, kappathr);
                dNh[jz][jM][0] = CLIGHT*PI*pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
                
                Mb = 2.0*PI*pow(rmax,2.0)*(C.dc(zl)-C.dc(zl-dz))*C.rhoM0;
                sigmab = interpolate(Mb, C.sigmalist);
                
                // see Baumann (5.129)
                dNh[jz][jM][1] = C.Dg(zl)*sigmab*C.halobias(zl, sigma);
                
                dNh[jz][jM][2] = rmax;
            }
        }
    }
    return dNh;
}

// number of filaments with kappa_CYL > kappa_thr
double NhfCYL(cosmology &C, double zs, double kappathr) {
    double Nh = 0.0;
    double zl, dz, M, dlnM, dndlnM, rmax;
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.FMFlist[jz][jM];
                rmax = rmaxfCYL(C, zs, zl, M, kappathr);
                Nh += CLIGHT*PI*pow((1.0+zl)*rmax,2.0)/C.Hz(zl)*dndlnM*dlnM*dz;
            }
        }
    }
    return Nh;
}


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                    PDF of amplifications                                                                       */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                       Correlated 1D environment field for the bias layer (bias_model = 1)                                       */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
//
// Replaces the legacy iid per-(jz,jM) log-normal count modulation with segment
// averages of ONE Gaussian field delta_1D(chi) along the LOS (ported from the
// emulator; design: docs/bias_field_design_note.md).
//
//   P_1D(kpar; Rperp) = (1/2pi) int dkperp kperp P0(sqrt(kpar^2+kperp^2))
//                       * W(k Rperp)^2            (KP91 eq. 3.8; window on k_perp
//   for the disk shape, on |k| for top-hat/Gaussian)
//   modes k_n = 2 pi n / L, L = 1.05 chi(z_s), n = 1..N_max = L/Rperp (>= 4);
//   per-mode field variance 2 P_1D(k_n)/L.
//
// The per-z-shell SEGMENT AVERAGES of the truncated mode sum form a Gaussian
// vector with covariance
//   Cov_ij = sum_n (2/L) P_1D(k_n) cos(k_n (c_i - c_j)) snc_i(k_n) snc_j(k_n),
//   snc_i(k) = sin(k L_i / 2)/(k L_i / 2),
// realized exactly by a Cholesky factor: dbar = chol * g, g ~ N(0,1)^n. P0(k) is
// the GROWTH-FREE linear power (cosmology::Pk0); growth enters per cell via
// b(M,z_l) Dg(z_l), exactly as in the legacy layer (Baumann 5.129 convention).

// Per-cell sub-threshold Campbell moments for the weak arm (bias_weak):
//   m[jz][jM] = int nbar kappa,  v[jz][jM] = int nbar kappa^2
// over the sub-threshold annuli, with the SAME integrand, log-annulus measure,
// stepping and floor as sigmakappaW so that sum_{jz,jM} v == sigmakappaW^2
// (asserted at build time; the only difference is per-cell bookkeeping).
static void weakMomentsNFW(cosmology &C, double zs, double kappathr, double eps_floor,
                           vector<vector<double> > &m, vector<vector<double> > &v) {
    m.assign(C.Nz, vector<double>(C.NM, 0.0));
    v.assign(C.Nz, vector<double>(C.NM, 0.0));
    double zl, dz, M, dlnM, dndlnM, r, kappar;
    double dlnr = 0.01;
    double Edlnr = exp(dlnr);
    for (int jz = 1; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        dz = zl - C.zlist[jz-1];
        if (zl < zs) {
            for (int jM = 1; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                dlnM = log(M) - log(C.Mlist[jM-1]);
                dndlnM = C.HMFlist[jz][jM][0];

                r = rmaxfNFW(C, zs, zl, M, kappathr);
                if (r == 0.0) {
                    r = 1.0e-6;
                }

                double SigmacW = Sigmacf(C, zs, zl);
                vector<double> NFWpW = interpolate2(zl, M, C.zlist, C.Mlist, C.NFWlist);
                double rsW = NFWpW[0];
                double kappa0W = kappa0NFW(rsW, NFWpW[1], SigmacW);

                kappar = kappathr;
                while (kappar > eps_floor*kappathr) {
                    kappar = kappagammaNFWeps(0.0, kappa0W, r/rsW, 0.0)[0];
                    double pref = CLIGHT*2.0*PI*pow((1.0+zl)*r,2.0)/C.Hz(zl)*dndlnM*dlnr*dlnM*dz;
                    m[jz][jM] += pref*kappar;
                    v[jz][jM] += pref*pow(kappar,2.0);
                    r = r*Edlnr;
                }
            }
        }
    }
}

// Smoothing window W~(x)^2 for the ISOTROPIC bias-field windows (bias_window
// 1 = spherical top-hat, 2 = Gaussian), x = |k| R. Window 0 (transverse disk)
// acts on k_perp alone and stays precomputed on the k_perp grid in build().
static inline double biasWindow2(double x, int window) {
    if (window == 2) {                                 // Gaussian, exp(-x^2/2)
        double W = exp(-0.5*x*x);
        return W*W;
    }
    // spherical top-hat, 3 (sin x - x cos x)/x^3; use the series below 1e-2 to
    // avoid cancellation.
    double W;
    if (x < 1.0e-2) {
        double x2 = x*x;
        W = 1.0 - x2/10.0*(1.0 - x2/28.0);
    } else {
        W = 3.0*(sin(x) - x*cos(x))/(x*x*x);
    }
    return W*W;
}

struct BiasField1D {
    int n = 0;                        // number of z-shells (jz = 1 .. n, zlist[jz] < zs)
    double Rperp = 0.0, L = 0.0;
    int window = 0;                   // cfg.bias_window (0 disk, 1 top-hat, 2 Gaussian)
    long Nmax = 0;                    // mode count L/Rperp (>= 4)
    std::vector<double> sig2;         // Cov_ii per shell (z=0 field, segment-averaged)
    std::vector<double> chol;         // lower-triangular Cholesky of Cov, row-major n*n

    // ---- weak (sub-threshold) arm tables, built only when bias_weak: per shell
    // i, on a uniform delta grid over +-DGRID_SIG sigma_i, log tables of
    // T_i(delta) = sum_M m_iM lambda_iM(delta) and V_i(delta) = sum_M v_iM
    // lambda_iM(delta); S_i = T_i - msum_i. Outside the grid delta is clamped.
    static constexpr int NGRID_W = 193;
    static constexpr double DGRID_SIG = 6.0;
    bool has_weak = false;
    std::vector<double> msum;         // per shell: sum_M m_iM
    std::vector<double> lnT, lnV;     // row-major n*NGRID_W

    // shell index for a given jz (cells at jz have chi in [dc(z_{jz-1}), dc(z_jz)])
    inline int shell(int jz) const { return (jz >= 1 && jz <= n) ? jz - 1 : -1; }

    // weak-arm lookup: conditional mean shift S and Campbell variance V of
    // shell i at field value delta (linear interp of the log tables)
    inline void weakSV(int i, double delta, double &S, double &V) const {
        double si = sqrt(sig2[i]);
        double half = DGRID_SIG*si;
        double x = std::min(std::max(delta, -half), half);
        double t = (x + half)/(2.0*half)*(NGRID_W - 1);
        int g = std::min(static_cast<int>(t), NGRID_W - 2);
        double f = t - g;
        const double *rT = &lnT[static_cast<size_t>(i)*NGRID_W];
        const double *rV = &lnV[static_cast<size_t>(i)*NGRID_W];
        S = exp((1.0 - f)*rT[g] + f*rT[g+1]) - msum[i];
        V = exp((1.0 - f)*rV[g] + f*rV[g+1]);
    }

    // Build the weak-arm tables. skappaW = sigmakappaW(...) with the SAME
    // kappathr/eps_floor — used for the sum_v == sigma_W^2 consistency gate.
    void buildWeak(cosmology &C, double zs, double kappathr, double eps_floor,
                   double skappaW) {
        if (n == 0) return;
        vector<vector<double> > mc, vc;
        weakMomentsNFW(C, zs, kappathr, eps_floor, mc, vc);
        double sv = 0.0;
        for (int jz = 0; jz < C.Nz; jz++)
            for (int jM = 0; jM < C.NM; jM++) sv += vc[jz][jM];
        if (fabs(sv - skappaW*skappaW) > 1.0e-9*skappaW*skappaW) {
            throw std::runtime_error("bias_weak: sum of per-cell v moments != sigma_W^2 "
                                     "(weakMomentsNFW drifted from sigmakappaW)");
        }
        msum.assign(n, 0.0);
        lnT.assign(static_cast<size_t>(n)*NGRID_W, 0.0);
        lnV.assign(static_cast<size_t>(n)*NGRID_W, 0.0);
        for (int i = 0; i < n; i++) {
            int jz = i + 1;
            double zl = C.zlist[jz];
            double Dgz = C.Dg(zl);
            double si = sqrt(sig2[i]);
            for (int g = 0; g < NGRID_W; g++) {
                double delta = (-DGRID_SIG + 2.0*DGRID_SIG*g/(NGRID_W - 1))*si;
                double T = 0.0, V = 0.0;
                for (int jM = 1; jM < C.NM; jM++) {
                    double mm = mc[jz][jM];
                    if (mm <= 0.0 && vc[jz][jM] <= 0.0) continue;
                    double a = Dgz*C.halobias(zl, C.sigmalist[jM][1]);
                    double lam = exp(a*delta - 0.5*a*a*sig2[i]);
                    T += mm*lam;
                    V += vc[jz][jM]*lam;
                }
                lnT[static_cast<size_t>(i)*NGRID_W + g] = log(std::max(T, 1.0e-300));
                lnV[static_cast<size_t>(i)*NGRID_W + g] = log(std::max(V, 1.0e-300));
            }
            for (int jM = 1; jM < C.NM; jM++) msum[i] += mc[jz][jM];
        }
        has_weak = true;
    }

    void build(cosmology &C, double zs, double Rp, int win = 0) {
        Rperp = Rp;
        window = win;
        // ---- shells
        std::vector<double> clo, chi_;
        for (int jz = 1; jz < C.Nz && C.zlist[jz] < zs; jz++) {
            clo.push_back(C.dc(C.zlist[jz-1]));
            chi_.push_back(C.dc(C.zlist[jz]));
        }
        n = static_cast<int>(clo.size());
        if (n == 0) return;

        L = 1.05*C.dc(zs);
        Nmax = std::max(4L, static_cast<long>(std::floor(L/Rperp)));
        if (Nmax > 5000000L) {
            throw std::invalid_argument("bias_Rperp too small: N_max = L/Rperp > 5e6 modes");
        }
        const double kmin = 2.0*PI/L;
        const double kmax = 2.0*PI*static_cast<double>(Nmax)/L;

        // ---- P_1D table (log-log interpolated; window 0 = disk via GSL J1 on
        // k_perp, precomputed here; windows 1/2 act on |k| and are evaluated
        // inside the k_par loop below).
        const int nkperp = 2048, nktab = 600;
        const double Rw = std::max(Rperp, 10.0);          // window floor, as in the prototype
        const double kperp_lo = 1.0e-9, kperp_hi = 60.0/Rw;
        const double dlnkp = log(kperp_hi/kperp_lo)/(nkperp - 1);
        std::vector<double> kperp(nkperp), W2(nkperp);
        for (int i = 0; i < nkperp; i++) {
            kperp[i] = kperp_lo*exp(dlnkp*i);
            double x = kperp[i]*Rperp;
            double Wd = (x < 1.0e-6) ? 1.0 : 2.0*gsl_sf_bessel_J1(x)/x;
            W2[i] = Wd*Wd;
        }
        std::vector<double> lktab(nktab), lPtab(nktab);
        const double lk0 = log(0.5*kmin), lk1 = log(kmax);
        for (int t = 0; t < nktab; t++) {
            double lk = lk0 + (lk1 - lk0)*t/(nktab - 1);
            double kpar = exp(lk);
            double s = 0.0, fprev = 0.0;
            for (int i = 0; i < nkperp; i++) {
                double kk = sqrt(kpar*kpar + kperp[i]*kperp[i]);
                double w2 = (window == 0) ? W2[i] : biasWindow2(kk*Rperp, window);
                double f = kperp[i]*kperp[i]*C.Pk0(kk)*w2;
                if (i > 0) s += 0.5*(f + fprev)*dlnkp;
                fprev = f;
            }
            lktab[t] = lk;
            lPtab[t] = log(std::max(s/(2.0*PI), 1.0e-300));
        }
        auto P1D = [&](double k) {
            double lk = log(k);
            if (lk <= lktab.front()) return exp(lPtab.front());
            if (lk >= lktab.back())  return exp(lPtab.back());
            int t = static_cast<int>((lk - lktab.front())/(lktab[1] - lktab[0]));
            t = std::min(t, nktab - 2);
            double f = (lk - lktab[t])/(lktab[t+1] - lktab[t]);
            return exp((1.0 - f)*lPtab[t] + f*lPtab[t+1]);
        };

        // ---- covariance: EXACT mode sum, k_q = 2 pi q / L, q = 1..Nmax.
        // Phases at the shell EDGES advance by a fixed rotation per mode
        // (uniform k grid); resynced every 4096 modes against drift.
        std::vector<double> a(n), b(n), Lh(n);
        for (int i = 0; i < n; i++) {
            a[i]  = clo[i];
            b[i]  = chi_[i];
            Lh[i] = chi_[i] - clo[i];
        }
        std::vector<double> Cov(static_cast<size_t>(n)*n, 0.0);
        std::vector<double> ca(n), sa(n), cb(n), sb(n);       // phases at edges
        std::vector<double> ra_c(n), ra_s(n), rb_c(n), rb_s(n); // per-mode rotations
        const double dk = 2.0*PI/L;
        for (int i = 0; i < n; i++) {
            ra_c[i] = cos(dk*a[i]); ra_s[i] = sin(dk*a[i]);
            rb_c[i] = cos(dk*b[i]); rb_s[i] = sin(dk*b[i]);
        }
        std::vector<double> cq(n), sq(n);
        for (long q = 1; q <= Nmax; q++) {
            double kq = dk*static_cast<double>(q);
            if (q == 1 || (q & 4095) == 0) {                  // init / resync
                for (int i = 0; i < n; i++) {
                    ca[i] = cos(kq*a[i]); sa[i] = sin(kq*a[i]);
                    cb[i] = cos(kq*b[i]); sb[i] = sin(kq*b[i]);
                }
            } else {                                          // rotate by dk
                for (int i = 0; i < n; i++) {
                    double c0 = ca[i], s0 = sa[i];
                    ca[i] = c0*ra_c[i] - s0*ra_s[i];
                    sa[i] = s0*ra_c[i] + c0*ra_s[i];
                    c0 = cb[i]; s0 = sb[i];
                    cb[i] = c0*rb_c[i] - s0*rb_s[i];
                    sb[i] = s0*rb_c[i] + c0*rb_s[i];
                }
            }
            double rw = sqrt(2.0*P1D(kq)/L);                  // per-mode field std
            for (int i = 0; i < n; i++) {
                double inv = 1.0/(kq*Lh[i]);
                cq[i] = rw*(sb[i] - sa[i])*inv;               // cos(k c) snc
                sq[i] = rw*(ca[i] - cb[i])*inv;               // sin(k c) snc
            }
            for (int i = 0; i < n; i++) {
                double *row = &Cov[static_cast<size_t>(i)*n];
                double ci = cq[i], si = sq[i];
                for (int jj = 0; jj <= i; jj++)
                    row[jj] += ci*cq[jj] + si*sq[jj];
            }
        }
        for (int i = 0; i < n; i++)
            for (int jj = i + 1; jj < n; jj++)
                Cov[static_cast<size_t>(i)*n + jj] = Cov[static_cast<size_t>(jj)*n + i];

        sig2.assign(n, 0.0);
        double trace = 0.0;
        for (int i = 0; i < n; i++) {
            sig2[i] = Cov[static_cast<size_t>(i)*n + i];
            trace += sig2[i];
        }

        // ---- Cholesky (PSD by construction; tiny relative jitter for roundoff)
        const double jitter = 1.0e-12*trace/n;
        for (int i = 0; i < n; i++) Cov[static_cast<size_t>(i)*n + i] += jitter;
        chol.assign(static_cast<size_t>(n)*n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int jj = 0; jj <= i; jj++) {
                double s = Cov[static_cast<size_t>(i)*n + jj];
                for (int kk = 0; kk < jj; kk++)
                    s -= chol[static_cast<size_t>(i)*n + kk]*chol[static_cast<size_t>(jj)*n + kk];
                if (i == jj) {
                    chol[static_cast<size_t>(i)*n + i] = sqrt(std::max(s, 0.0));
                } else {
                    double d = chol[static_cast<size_t>(jj)*n + jj];
                    chol[static_cast<size_t>(i)*n + jj] = (d > 0.0) ? s/d : 0.0;
                }
            }
        }
    }
};


// find threshold kappa
double findkappathr(int N, function<double(double)> Nf) {
    double kappa1 = 1.0e-12, kappa2 = 1.0;
    double kappathr = pow(10.0, (log10(kappa1) + log10(kappa2))/2.0);
    while (log10(kappa2) - log10(kappa1) > 0.01) {
        if (Nf(kappathr) > N) {
            kappa1 = kappathr;
        } else {
            kappa2 = kappathr;
        }
        kappathr = pow(10.0, (log10(kappa1) + log10(kappa2))/2.0);
    }
    return kappathr;
}


// probability distribution of lnmu, {lnmu, dP/dlnmu}
vector<vector<double> > lensing::Plnmuf(cosmology &C, double zs, rgen &mt, int fil, int bias, int ell, int write) {

    // ---- clustering-bias configuration guards (see lensing.h)
    if (bias_model != 0 && bias_model != 1) {
        throw std::invalid_argument("bias_model must be 0 (legacy iid cell bias) or 1 (correlated 1D field)");
    }
    if (bias_model == 1 && bias_Rperp <= 0.0) {
        throw std::invalid_argument("bias_model = 1 requires bias_Rperp > 0 (comoving kpc)");
    }
    if (bias_window < 0 || bias_window > 2) {
        throw std::invalid_argument("bias_window must be 0 (transverse disk), 1 (spherical top-hat) or 2 (Gaussian)");
    }
    if (bias_window != 0 && bias_model != 1) {
        throw std::invalid_argument("bias_window != 0 requires bias_model = 1 (the window shapes the correlated field's power spectrum)");
    }
    if (bias_weak && bias_model != 1) {
        throw std::invalid_argument("bias_weak requires bias_model = 1 (the correlated field supplies the conditioning); it is a no-op when bias = 0");
    }
    if (fil_bias && bias_model != 1) {
        throw std::invalid_argument("fil_bias requires bias_model = 1 (the correlated field carries the filament modulation); it is a no-op in the legacy iid layer");
    }

    // ---- subhalo configuration guards (see lensing.h)
    if (subhalo && subhalo_model != 1 && subhalo_model != 3
        && subhalo_model != 4 && subhalo_model != 5) {
        throw std::invalid_argument("subhalo_model must be 1 (reduced-host resolved-only), 3 (+ Gaussian unresolved term), 4 (brute) or 5 (kappa-thresholded brute)");
    }
    if (subhalo && subhalo_model == 3 && subhalo_brute) {
        throw std::invalid_argument("subhalo_brute is incompatible with subhalo_model 3 (use model 1 + subhalo_brute + subhalo_carve for a brute reference)");
    }
    // model 1 reduces the host by the resolved fraction f_s,res(r), which here is
    // only formed along the carve path. The fsb tables it would otherwise need
    // are built for model 3 only, so without the carve the host would be left
    // unreduced and the substructure mass counted twice
    if (subhalo && subhalo_model == 1 && !subhalo_carve) {
        throw std::invalid_argument("subhalo_model 1 requires subhalo_carve = true in this port (the non-carve host reduction by f_s,res(r) is not implemented here; use model 3 for the deterministic reduction)");
    }
    if (subhalo && subhalo_model == 4 && !subhalo_carve) {
        throw std::invalid_argument("subhalo_model 4 requires subhalo_carve = true (the carve is intrinsic to the model)");
    }
    if (subhalo && subhalo_model == 5 && !subhalo_carve) {
        throw std::invalid_argument("subhalo_model 5 requires subhalo_carve = true (the carve is intrinsic to the model)");
    }
    if (subhalo && subhalo_model == 5 && subhalo_brute) {
        throw std::invalid_argument("subhalo_brute is meaningless for subhalo_model 5 (the model is brute by construction, thresholded on kappa)");
    }
    if (subhalo && subhalo_model == 5 && subhalo_kappathr <= 0.0
        && subhalo_kappathr_factor <= 0.0) {
        throw std::invalid_argument("subhalo_model 5 needs subhalo_kappathr > 0 or subhalo_kappathr_factor > 0");
    }
    if (subhalo && subhalo_virial
        && !(subhalo_model == 4 || subhalo_model == 5)) {
        throw std::invalid_argument("subhalo_virial requires subhalo_model 4 or 5 (models 1/3 reduce the host with M_200-referred incomplete-Gamma/Wsub tables)");
    }

    // fix threshold kappa
    function<double(double)> NfNFW = [&C, zs](double kappa) {
        return NhfNFW(C, zs, kappa);
    };
    double kappathrH = findkappathr(Nhalos, NfNFW);
    
    if (write > 0) {
        cout << kappathrH << endl;
    }
        
    // distribution of kappa_NFW < kappa_thr
    double skappaW = sigmakappaW(C, zs, kappathrH);
    normal_distribution<double> PkappaW(0.0, skappaW);
    if (skappaW < 0.0) {
        cout << "Error: negative standard deviation." << endl;
    }
    
    // bias_weak (+ field + bias on): the weak background is drawn CONDITIONALLY on
    // the field, which is only available after the field block below — the initial
    // fill is skipped (0.0) here and done there instead. All other paths keep the
    // legacy unconditional fill (bit-identical stream for bias_model = 0).
    const bool weak_conditional = bias_weak && bias_model == 1 && bias != 0;

    vector<double> kappalist(Nreal, 0.0);
    vector<double> gamma1list(Nreal, 0.0);
    vector<double> gamma2list(Nreal, 0.0);
    for (int j = 0; j < Nreal; j++) {
        kappalist[j] = weak_conditional ? 0.0 : PkappaW(mt);
    }

    vector<vector<vector<double> > > dNH = deltaNhfNFW(C, zs, kappathrH);
    vector<vector<vector<double> > > dNF = deltaNhfCYL(C, zs, kappathrH);
    if (subhalo) {
        S.m_floor = subhalo_m_floor;
        S.psi_min_fixed = psi_min_fixed;
        S.virial = subhalo_virial; // set before precompute
        // model 3 also builds the Wsub tables, which need the host threshold for
        // the encounter disc radius. Model 5 reuses r_thr as the clump reach
        // D(m), so there it must be built at kappathr_sub instead
        double clump_kappathr = subhalo_factor*kappathrH;
        if (subhalo_model == 5) {
            clump_kappathr = (subhalo_kappathr > 0.0)
                               ? subhalo_kappathr
                               : subhalo_kappathr_factor*kappathrH;
        }
        S.precompute(C, zs, clump_kappathr,
                     (subhalo_model == 3) ? kappathrH : 0.0,
                     subhalo_model == 5);
    }
    if (write > 0) {
        writeToFile(C.zlist, C.Mlist, dNH, C.outdir/"dNH.dat");
        writeToFile(C.zlist, C.Mlist, dNF, C.outdir/"dNF.dat");
    }
    
    array<double,2> kappagamma;
    normal_distribution<double> pG(0.0, 1.0);

    // bias_model = 1: draw the correlated per-shell environment field for ALL
    // realizations up front (the sampling loops below are cell-major, so shell jz
    // needs every realization's field value when its cells are processed). This is
    // the ONLY model-1 RNG consumption outside the shared path; model 0 draws
    // nothing here and is bit-identical to the legacy stream. bias == 0 skips the
    // build (lambda is forced to 1, so the field would be dead weight and its draws
    // would needlessly shift the halo stream vs bias_model = 0).
    // Memory: Nreal * n_shells floats.
    BiasField1D bfield;
    std::vector<float> bfvals;
    if (bias_model == 1 && bias != 0) {
        bfield.build(C, zs, bias_Rperp, bias_window);
        if (bfield.n > 0) {
            const int nsh = bfield.n;
            bfvals.resize(static_cast<size_t>(Nreal)*nsh);
            std::vector<double> g(nsh);
            for (int j = 0; j < Nreal; j++) {
                for (int i = 0; i < nsh; i++) g[i] = pG(mt);
                for (int i = 0; i < nsh; i++) {
                    double s = 0.0;
                    const double *row = &bfield.chol[static_cast<size_t>(i)*nsh];
                    for (int kk = 0; kk <= i; kk++) s += row[kk]*g[kk];
                    bfvals[static_cast<size_t>(j)*nsh + i] = static_cast<float>(s);
                }
            }
        }
    }

    // bias_weak: conditional weak background, drawn from the SAME realized field
    // values as the count modulation (Cox split of Campbell's theorem). Consumes
    // one normal per realization, at this fixed stream point. The eps floor 0.001
    // matches sigmakappaW's hardcoded 0.001*kappathr, so the sum_v == sigma_W^2
    // gate in buildWeak holds. Degenerate n == 0 (z_s below the first grid shell)
    // falls back to the legacy unconditional draw.
    if (weak_conditional) {
        if (bfield.n > 0) {
            bfield.buildWeak(C, zs, kappathrH, 0.001, skappaW);
            const int nsh = bfield.n;
            for (int j = 0; j < Nreal; j++) {
                double sumS = 0.0, sumV = 0.0, Sw, Vw;
                const float *fj = &bfvals[static_cast<size_t>(j)*nsh];
                for (int i = 0; i < nsh; i++) {
                    bfield.weakSV(i, static_cast<double>(fj[i]), Sw, Vw);
                    sumS += Sw;
                    sumV += Vw;
                }
                kappalist[j] = sumS + sqrt(std::max(sumV, 0.0))*pG(mt);
            }
        } else {
            for (int j = 0; j < Nreal; j++) {
                kappalist[j] = PkappaW(mt);
            }
        }
    }

    poisson_distribution<int> PN;
    double zl, M, rmaxH, rmaxF, r, phi, phiH, phiF, epsilon = 0.0, barNH, barNF, sigma, deltab, lambda, meankappa = 0.0;
    int NH, NF;
    int NtotH = 0, NtotF = 0;;

    // whether the mass conserving carve path applies
    const bool carve_path = subhalo && subhalo_carve &&
        (subhalo_model == 3 || subhalo_model == 4 || subhalo_model == 5 ||
         (subhalo_brute && subhalo_model == 1));

    // one host encounter with substructure, shared by the two structurally
    // identical branches of the loop below, the small lambda shortcut and the
    // Poisson branch, so that they cannot diverge
    auto subhaloEncounter = [&](int j, int jz, int jM, double zl_, double M_,
                                double Sigmac_, double r_, double phi_, double phiH_) {
        const double kappabefore = kappalist[j];
        array<double,2> kg;
        if (carve_path) {
            // clumps first, accumulating their realized mass. For model 5 that
            // is the retained mass and the rest stays in the smooth host, which
            // is what makes the threshold mass conserving. The host build draws
            // no random numbers, so this order leaves the stream unchanged
            double Msum = 0.0;
            if (subhalo_model == 5) {
                S.addClumpsRestricted(C, jz, jM, M_, Sigmac_, r_, phi_, mt,
                                      kappalist[j], gamma1list[j], gamma2list[j], &Msum);
            } else {
                S.addClumps(C, jz, jM, zl_, M_, Sigmac_, r_, phi_, mt,
                            kappalist[j], gamma1list[j], gamma2list[j],
                            subhalo_model, subhalo_brute, &Msum);
            }

            // mean unresolved mass at this ray, model 3 only
            const double M_u = (subhalo_model == 3)
                                 ? S.unresolvedMass(C, jz, jM, r_, M_) : 0.0;

            // carved host at M - sum_i m_i - M_u, with a negative mass guard. In
            // virial mode Msum is drawn against the M_vir budget while the host
            // is parameterized by its M200 grid mass, so convert: the host then
            // keeps the same fractional mass 1-f_s in both apertures, and the
            // remainder is substructure in the r200 to r_vir shell, which was
            // never part of the M200 budget
            const double vr = S.virialRatio(jz, jM);
            double Mhost = M_ - Msum/vr - M_u;
            if (Mhost < C.Mmin) Mhost = C.Mmin;

            vector<double> NFWp = interpolate2(zl_, Mhost, C.zlist, C.Mlist, C.NFWlist);
            double rs = NFWp[0];
            double kappa0 = kappa0NFW(rs, NFWp[1], Sigmac_);
            double eps = (ell > 0) ? epsilonNFW(C, zl_, Mhost) : 0.0;
            kg = kappagammaNFWeps(eps, kappa0, r_/rs, phiH_);
            kappalist[j] += kg[0];
            gamma1list[j] += cos(phi_)*kg[1];
            gamma2list[j] += sin(phi_)*kg[1];

            // unresolved clump term, model 3 only, at the same position in the
            // random stream as on the non-carve path below
            if (subhalo_model == 3) {
                double muU, sU;
                S.wsubTerm(jz, jM, r_, muU, sU);
                kappalist[j] += muU + sU*pG(mt);
            }
        } else {
            // deterministic reduction, the host at (1-f_s,b)M, reduced by the
            // mean bound fraction rather than the realized clump mass
            double Msm = max(C.Mlist[0], (1.0 - S.fsb[jz][jM])*M_);
            vector<double> NFWp = interpolate2(zl_, Msm, C.zlist, C.Mlist, C.NFWlist);
            double rs = NFWp[0];
            double kappa0 = kappa0NFW(rs, NFWp[1], Sigmac_);
            double eps = (ell > 0) ? epsilonNFW(C, zl_, Msm) : 0.0;
            kg = kappagammaNFWeps(eps, kappa0, r_/rs, phiH_);
            kappalist[j] += kg[0];
            gamma1list[j] += cos(phi_)*kg[1];
            gamma2list[j] += sin(phi_)*kg[1];
            S.addClumps(C, jz, jM, zl_, M_, Sigmac_, r_, phi_, mt,
                        kappalist[j], gamma1list[j], gamma2list[j],
                        subhalo_model, subhalo_brute, nullptr);
            if (subhalo_model == 3) {
                double muW, sigmaW;
                S.wsubTerm(jz, jM, r_, muW, sigmaW);
                kappalist[j] += muW + sigmaW*pG(mt);
            }
        }
        meankappa += kappalist[j] - kappabefore;
    };
    for (int jz = 0; jz < C.Nz; jz++) {
        zl = C.zlist[jz];
        if (zl < zs) {
            for (int jM = 0; jM < C.NM; jM++) {
                M = C.Mlist[jM];
                
                // generate realizations
                barNH = dNH[jz][jM][0];
                sigma = dNH[jz][jM][1];
                rmaxH = dNH[jz][jM][2];
                
                barNF = dNF[jz][jM][0];
                rmaxF = dNF[jz][jM][2];

                double Sigmac = (barNH > 0.0 || barNF > 0.0) ? Sigmacf(C, zs, zl) : 0.0;
                double rsH = 1.0, kappa0H = 0.0;
                if (barNH > 0.0) {
                    vector<double> NFWp = interpolate2(zl, M, C.zlist, C.Mlist, C.NFWlist);
                    rsH = NFWp[0];
                    kappa0H = kappa0NFW(rsH, NFWp[1], Sigmac);
                    if (ell > 0) epsilon = epsilonNFW(C, zl, M);
                }
                double rsF = 0.0, LF = 0.0, kappa0baseF = 0.0;
                if (fil > 0 && barNF > 0.0) {
                    rsF = 1000.0*pow(M/1.0e14, 1.0/3.0);
                    LF = 20000.0*pow(M/1.0e14, 1.0/3.0);
                    kappa0baseF = rsF * 14.4*C.rhoc*LF / Sigmac;
                }

                // bias_model = 1: this cell's shell index and field amplitudes.
                // bDg = b(M,z_l) Dg(z_l) (growth-free field; the sigmab variance is
                // carried by the field's own sig2), bcomp = 1/2 bDg^2 sig2 so that
                // <lambda> = 1 exactly. bDgF/bcompF use the filament bias filbias
                // (cfg fil_bias) so filaments ride b_fil instead of the halo bias;
                // they stay 0 (=> lambdaF falls back to lambda) otherwise.
                int bsh = -1;
                double bDg = 0.0, bcomp = 0.0, bDgF = 0.0, bcompF = 0.0;
                if (bias_model == 1) {
                    bsh = bfield.shell(jz);
                    if (bsh >= 0) {
                        bDg = C.Dg(zl)*C.halobias(zl, C.sigmalist[jM][1]);
                        bcomp = 0.5*bDg*bDg*bfield.sig2[bsh];
                        if (fil_bias) {
                            bDgF = C.Dg(zl)*C.filbias(zl, C.sigmalist[jM][1]);
                            bcompF = 0.5*bDgF*bDgF*bfield.sig2[bsh];
                        }
                    }
                }

                for (int j = 0; j < Nreal; j++) {

                    // bias: correlated field (bias_model = 1) or legacy iid draw
                    if (bias_model == 1) {
                        lambda = (bsh >= 0)
                            ? exp(bDg*bfvals[static_cast<size_t>(j)*bfield.n + bsh] - bcomp)
                            : 1.0;
                    } else {
                        deltab = sigma*pG(mt);
                        lambda = exp(deltab - pow(sigma,2.0)/2.0); // log-normal
                    }

                    if (bias == 0) {
                        lambda = 1.0;
                    }

                    // Filament count modulation. Same realized field, but the
                    // filament amplitude uses filbias (fil_bias, bias_model = 1).
                    // Defaults to the halo lambda, so fil_bias = 0 (and the whole
                    // legacy layer) is bitwise unchanged; no new RNG draw.
                    double lambdaF = lambda;
                    if (fil_bias && bias != 0 && bias_model == 1 && bsh >= 0) {
                        lambdaF = exp(bDgF*bfvals[static_cast<size_t>(j)*bfield.n + bsh] - bcompF);
                    }

                    // generate halos
                    if (lambda*barNH < 0.2) { // if lambda is small, compare to a random number U(0,1) (faster)
                        if (lambda*barNH > randomreal(0.0, 1.0, mt)) {
                            r = sqrt(randomreal(0.0,1.0,mt))*rmaxH; // distance from the line-of-sight
                            phi = randomreal(0.0,2*PI,mt); // polar angle of r vector
                            phiH = randomreal(0.0,2*PI,mt); // orientation of the halo ellipticity

                            if (subhalo) {
                                subhaloEncounter(j, jz, jM, zl, M, Sigmac, r, phi, phiH);
                            } else {
                                kappagamma = kappagammaNFWeps(epsilon, kappa0H, r/rsH, phiH);

                                kappalist[j] += kappagamma[0];
                                gamma1list[j] += cos(phi)*kappagamma[1];
                                gamma2list[j] += sin(phi)*kappagamma[1];

                                meankappa += kappagamma[0];
                            }
                            
                            NtotH++;
                        }
                    } else { // for larger lambda, generate number of halos from Poisson distribution (slower)
                        PN = poisson_distribution<int>(lambda*barNH);
                        NH = PN(mt);
                        if (NH > 0) {
                            for (int jH = 0; jH < NH; jH++) {
                                r = sqrt(randomreal(0.0,1.0,mt))*rmaxH; // distance from the line-of-sight
                                phi = randomreal(0.0,2*PI,mt); // polar angle of r vector
                                phiH = randomreal(0.0,2*PI,mt); // orientation of the halo ellipticity

                                if (subhalo) {
                                    subhaloEncounter(j, jz, jM, zl, M, Sigmac, r, phi, phiH);
                                } else {
                                    kappagamma = kappagammaNFWeps(epsilon, kappa0H, r/rsH, phiH);
                                    
                                    kappalist[j] += kappagamma[0];
                                    gamma1list[j] += cos(phi)*kappagamma[1];
                                    gamma2list[j] += sin(phi)*kappagamma[1];
                                    
                                    meankappa += kappagamma[0];
                                }
                                
                                NtotH++;
                            }
                        }
                    }
                    
                    if (fil > 0) {
                        // generate filaments
                        if (lambdaF*barNF < 0.2) { // if lambda is small, compare to a random number U(0,1) (faster)
                            if (lambdaF*barNF > randomreal(0.0, 1.0, mt)) {
                                r = sqrt(randomreal(0.0,1.0,mt))*rmaxF; // distance from the line-of-sight
                                phi = randomreal(0.0,2*PI,mt); // polar angle of r vector
                                phiF = randomreal(0.0,2*PI,mt); // orientation of the filament
                                double kappa0F = kappa0baseF / max(2.0*rsF, abs(cos(phiF))*LF);
                                kappagamma = {kappaCYL2(r, rsF, kappa0F), gammaCYL2(r, rsF, kappa0F)};
                                
                                kappalist[j] += kappagamma[0];
                                gamma1list[j] += cos(phi)*kappagamma[1];
                                gamma2list[j] += sin(phi)*kappagamma[1];
                                
                                meankappa += kappagamma[0];
                                
                                NtotF++;
                            }
                        } else { // for larger lambda, generate number of halos from Poisson distribution (slower)
                            PN = poisson_distribution<int>(lambdaF*barNF);
                            NF = PN(mt);
                            if (NF > 0) {
                                for (int jF = 0; jF < NF; jF++) {
                                    r = sqrt(randomreal(0.0,1.0,mt))*rmaxF; // distance from the line-of-sight
                                    phi = randomreal(0.0,2*PI,mt); // polar angle of r vector
                                    phiF = randomreal(0.0,2*PI,mt); // orientation of the filament
                                    double kappa0F = kappa0baseF / max(2.0*rsF, abs(cos(phiF))*LF);
                                    kappagamma = {kappaCYL2(r, rsF, kappa0F), gammaCYL2(r, rsF, kappa0F)};
                                    
                                    kappalist[j] += kappagamma[0];
                                    gamma1list[j] += cos(phi)*kappagamma[1];
                                    gamma2list[j] += sin(phi)*kappagamma[1];
                                    
                                    meankappa += kappagamma[0];
                                    
                                    NtotF++;
                                }
                            }
                        }
                    }
                    
                }
            }
        }
    }
    // anchor that enforces <kappa> = 0, see lensing.h. Modes 1 and 2 read the
    // per-ray totals, which include the realized weak background, whereas the
    // running sum does not, so they differ at order sigmaW/sqrt(Nreal)
    if (kappa_anchor == 2) {
        meankappa = kappa_anchor_value;
    } else if (kappa_anchor == 1) {
        // drop the rays above the cut so that no single one shifts the batch
        double sum = 0.0;
        long nk = 0;
        for (int j = 0; j < Nreal; j++) {
            if (kappalist[j] <= kappa_anchor_cut) { sum += kappalist[j]; nk++; }
        }
        if (nk > 0) {
            meankappa = sum/(1.0*nk);
        } else { // every ray above the cut, fall back to the running sum
            meankappa = meankappa/(1.0*Nreal);
        }
    } else {
        meankappa = meankappa/(1.0*Nreal);
    }

    //cout << NtotH/(1.0*Nreal) << "   " << NtotF/(1.0*Nreal) << endl;
    
    // compute mu
    vector<double> lnmulist, lnmuAlist, ln1pkappalist, lngammalist;
    double kappaj, gamma1j, gamma2j, gammaj, muj;
    for (int j = 0; j < Nreal; j++) {
        kappaj = kappalist[j] - meankappa;
        gammaj = sqrt(pow(gamma1list[j], 2.0) + pow(gamma2list[j], 2.0));
        muj = 1.0/(pow(1.0-kappaj, 2.0) - pow(gammaj, 2.0));
        
        if (muj > 0.0 && gammaj > 0.0) {
            lnmulist.push_back(log(muj));
            lnmuAlist.push_back(log(1.0/pow(1.0-kappaj, 2.0)));
            ln1pkappalist.push_back(log(1.0+kappaj));
            lngammalist.push_back(log(gammaj));
        }
    }
    
    // binning
    vector<vector<double> > Plnmu = binSample(lnmulist, Nbins);
    if (write > 0) {
        vector<vector<double> > PlnmuA = binSample(lnmuAlist, Nbins);
        vector<vector<double> > Pln1pkappa = binSample(ln1pkappalist, Nbins);
        vector<vector<double> > Plngamma = binSample(lngammalist, Nbins);
            
        writeToFile(Plnmu, C.outdir/("Plnmu_z=" + to_string_prec(zs,1) + ".dat"));
        writeToFile(PlnmuA, C.outdir/("PlnmuA_z=" + to_string_prec(zs,1) + ".dat"));
        writeToFile(Pln1pkappa, C.outdir/("Pln1pkappa_z=" + to_string_prec(zs,1) + ".dat"));
        writeToFile(Plngamma, C.outdir/("Plngamma_z=" + to_string_prec(zs,1) + ".dat"));
    }
    
    int jmin = 0, jmax = Plnmu.size()-1;
    for (int j = 0; j < Plnmu.size()-4; j++) {
        if (Plnmu[j][1] > 0 && Plnmu[j+1][1] > 0 && Plnmu[j+2][1] > 0 && Plnmu[j+3][1] > 0 && Plnmu[j+4][1] > 0) {
            jmin = j;
            j = Plnmu.size();
        }
    }
    for (int j = max(4,jmin); j < Plnmu.size(); j++) {
        if (Plnmu[j-4][1] == 0 && Plnmu[j-3][1] == 0 && Plnmu[j-2][1] == 0 && Plnmu[j-1][1] == 0 && Plnmu[j][1] == 0) {
            jmax = j;
            j = Plnmu.size();
        }
    }
    //cout << jmin << "   " << jmax << endl;
    Plnmu.erase(Plnmu.begin() + jmax, Plnmu.end());
    Plnmu.erase(Plnmu.begin(), Plnmu.begin() + jmin);
    
    // convert from image plane to source plane (P_S ~ P_I/mu) and normalize
    double dlnmu = Plnmu[1][0] - Plnmu[0][0];
    double norm = 0.0;
    for (int j = 0; j < Plnmu.size(); j++) {
        Plnmu[j][1] *= exp(-Plnmu[j][0]);
        norm += Plnmu[j][1]*dlnmu;
    }
    for (int j = 0; j < Plnmu.size(); j++) {
        Plnmu[j][1] *= 1.0/norm;
    }
    return Plnmu;
}


/* ---------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                      Likelihood analysis of the Hubble diagram                                                                 */
/* ---------------------------------------------------------------------------------------------------------------------------------------------- */

// loglikelihood of the Hubble digram data
double lensing::loglikelihood(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &par, int lens, int dm, rgen &mt) {

    // initialize cosmology
    C.OmegaM = par[0];
    C.sigma8 = par[1];
    C.h = par[2];
    C.initialize(dm, pow(10.0, par[3]));
    
    double z, DL0, DL, sigmaDL, Y, dY, Pdet;
    
    // compute loglikelihood
    double logL = 0.0;
    if (lens == 0) { // model without lensing
        for (int j = 0; j < data.size(); j++) {
            z = data[j][0];
            DL0 = C.DL(z);
            DL = data[j][1];
            sigmaDL = data[j][2];
            
            // compute P_det
            dY = sigmaDL/10.0;
            Y = min(DL0, DLthr) - 3.0*sigmaDL;
            Pdet = 0.0;
            while (Y <= min(DLthr, DL0 + 3.0*sigmaDL)) {
                Pdet += dY*NPDF(Y, DL0, sigmaDL);
                Y += dY;
            }
            logL += logNPDF(DL, DL0, sigmaDL) - log(Pdet);
        }
    }
    else { // model with lensing
        
        // compute the lensing distributions
        vector<vector<vector<double> > > Plnmuz(C.Zlist.size());
        for (int jz = 0; jz < C.Zlist.size(); jz++) {
            z = C.Zlist[jz];
            Plnmuz[jz] = Plnmuf(C, z, mt, 1, 1, 1, 0);
        }
        
        int jz;
        vector<vector<double> > Plnmu;
        double L, dlnmu;
        for (int j = 0; j < data.size(); j++) {
            z = data[j][0];
            
            DL = data[j][1];
            sigmaDL = data[j][2];
            
            // find the closest z at which P(lnmu) is computed
            jz = lower_bound(C.Zlist.begin(), C.Zlist.end(), z) - C.Zlist.begin();
            if ((jz > 0 && C.Zlist[jz]-z > z-C.Zlist[jz-1]) || jz >= C.Zlist.size()) {
                jz--;
            }
            Plnmu = Plnmuz[jz];
            dlnmu = Plnmu[1][0] - Plnmu[0][0];
            
            // integrate over lnmu
            dY = sigmaDL/10.0;
            Pdet = 0.0, L = 0.0;
            for (int i = 0; i < Plnmu.size(); i++) {
                DL0 = C.DL(z)/exp(Plnmu[i][0]/2.0);
                
                L += dlnmu*Plnmu[i][1]*NPDF(DL, DL0, sigmaDL);
                
                // integrate over Y
                Y = min(DL0, DLthr) - 3.0*sigmaDL;
                while (Y <= min(DLthr, DL0 + 3.0*sigmaDL)) {
                    Pdet += dlnmu*Plnmu[i][1]*dY*NPDF(Y, DL0, sigmaDL);
                    Y += dY;
                }
            }
            
            logL += log(L/Pdet);
        }
    }
    
    return logL;
}

// MCMC inference of the Hubble diagram data
void lensing::Hubble_diagram_fit(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int N, int Nburnin, int lens, int dm, rgen &mt, fs::path filename) {
    
    // loglikelihood
    function<double(vector<double>&)> logpdf = [this, &C, DLthr, &data, lens, dm, &mt](vector<double> &par) {
        return loglikelihood(C, DLthr, data, par, lens, dm, mt);
    };
    
    // no cut
    function<double(vector<double>&)> cut = [](vector<double> &par) {
        return 1.0;
    };
    
    MCMC_sampling(N, Nburnin, logpdf, initial, steps, priors, cut, mt, 1, 0, filename);
}
