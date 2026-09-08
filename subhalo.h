class cosmology;

class Subhalo {

public:
    // psi_max = 1 matches the range over which gden normalizes the SHMF; the
    // exp(-beta psi^omega) cutoff supplies the suppression near psi -> 1
    double alpha = -0.82, beta = 50.0, omega = 4.0, psi_res = 1.0e-4, psi_max = 1.0;
    double m_floor = 1.0e7;
    // replaces the absolute floor psi_min = m_floor/M with a fixed psi_min
    // everywhere it is computed, e.g. psi_res to avoid extrapolating the SHMF
    // below its own calibration point. <= 0 keeps m_floor/M
    double psi_min_fixed = -1.0;
    int Nu = 128, NyW = 48;

    // Jiang & van den Bosch 2014 define haloes and subhaloes inside their virial
    // radii, so f_s and psi = m/M are M_vir quantities and the population reaches
    // r_vir, and Green et al. 2021 normalize the radial bias at r_vir. The grid
    // mass here is M200c and the profile was sampled only to r200, which leaves
    // 4-20% too much substructure inside r200. virial refers psi to
    // M_vir = M mu(c_vir)/mu(c200) and samples out to x = eta = r_vir/r200.
    // Gated to subhalo_model 4 and 5, which throw otherwise
    bool virial = false;
    vector<vector<double> > xmaxh; // clump extent in r200 units: 1, or eta
    vector<vector<double> > Mpsih; // mass scale for psi = m/Mpsi: M, or M_vir
    vector<vector<double> > Mgrid_; // grid mass M per bin, for virialRatio

    // M_vir/M200 for the host bin, 1 in legacy mode. lensing.cpp converts the
    // realized clump mass to the M200 scale before carving, so that the smooth
    // host keeps the same fractional mass 1-f_s in both apertures
    double virialRatio(int jz, int jM) const {
        if (!virial || Mpsih.empty()) return 1.0;
        double Mv = Mpsih[jz][jM];
        return (Mv > 0.0) ? Mv/Mgrid_[jz][jM] : 1.0;
    }

    vector<vector<double> > gnorm;
    vector<vector<double> > r200h;
    vector<vector<double> > chost;
    vector<vector<double> > r_thr;
    vector<vector<vector<double> > > invRad;
    vector<vector<vector<double> > > muW;
    vector<vector<vector<double> > > sW;
    vector<vector<array<double,2> > > lyW;
    vector<vector<double> > fsb;
    double log_Mmin = 0.0;
    double inv_dlogM = 0.0;

    // restricted intensity sampling, subhalo_model 5. The population is the same
    // as model 4, every subhalo down to psi_min = m_floor/M with no unresolved
    // stand-in, but only clumps whose kappa at the ray clears kappathr_sub are
    // rendered, and that restricted Poisson intensity is sampled directly so the
    // rejects are never instantiated. Drawing and rejecting would save nothing,
    // since the test needs both m and d.
    //
    // The retention region is the disc d <= D(m) around the ray, with D the clump
    // reach r_thr built at kappathr_sub. Writing the projected clump number
    // density as Sigma_n(R), normalized so that int Sigma_n dA = 1, the retained
    // intensity per clump mass is int_disc Sigma_n dA. Thinning against the
    // global envelope Smax = max_R Sigma_n, propose with probability
    // min(1, pi D^2 Smax), then either place the clump uniformly in the disc and
    // accept with Sigma_n(R)/Smax, or, once that probability caps at 1, draw from
    // the full radial profile and accept if d <= D. Both branches are exact and
    // the envelope only costs extra proposals. Since Smax is global the proposal
    // intensity does not depend on the ray, so it precomputes per bin.
    // Dropped clumps keep their mass in the smooth host through the carve, so
    // mass_out accumulates retained clump mass only.
    //
    // This is the rule the host haloes already obey, applied to subhaloes, not
    // the weak and unresolved split of model 3.
    int NRp = 128; // projected profile radial bins
    int Nq = 96; // log-psi proposal grid points
    vector<vector<vector<double> > > p2p; // normalized projected profile
    vector<vector<double> > sig2dmax; // envelope max of Sigma_n
    vector<vector<vector<double> > > propCum; // cumulative proposal intensity
    vector<vector<array<double,2> > > lpq; // {log psi_lo, dlog psi}
    bool restricted_built = false;

    double sigma2Dclump(int jz, int jM, double s) const;
    double clumpReach(cosmology &C, int jz, double m) const;
    void buildRestrictedBin(cosmology &C, int jz, int jM);

    // kappathr_host > 0 also builds the Wsub tables; build_restricted also
    // builds the model 5 tables, and then kappathr must be kappathr_sub since
    // r_thr doubles as the clump reach D(m)
    void precompute(cosmology &C, double zs, double kappathr, double kappathr_host,
                    bool build_restricted = false);
    void buildWsubBin(cosmology &C, double zs, int jz, int jM, double kappathr_host);
    void wsubTerm(int jz, int jM, double r, double &mu, double &sigma);
    double unresolvedMass(cosmology &C, int jz, int jM, double r, double M) const;

    // mass_out, if given, accumulates the realized clump mass sum_i m_i that the
    // mass conserving carve removes from the host
    int addClumps(cosmology &C, int jz, int jM, double zl, double M, double Sigmac, double r, double phi, rgen &mt, double &kappa, double &gamma1, double &gamma2,
                  int subhalo_model = 3, bool subhalo_brute = false, double *mass_out = nullptr);

    // as addClumps, but draws only the clumps that clear kappathr_sub, and
    // returns the number rendered rather than proposed
    int addClumpsRestricted(cosmology &C, int jz, int jM, double M, double Sigmac,
                            double r, double phi, rgen &mt,
                            double &kappa, double &gamma1, double &gamma2,
                            double *mass_out = nullptr);

};
