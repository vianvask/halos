
#include "subhalo.h"

class lensing {
    
public:
    int Nreal; // realizations
    int Nhalos; // number of halos in each realization
    int Nbins; // P(lnmu) bins
    // ---- Subhalo substructure. The production settings are a coupled set, so a
    // single flag override can throw: subhalo_model = 3 needs subhalo_virial off,
    // and models 4 and 5 need subhalo_carve on.
    bool subhalo = true; // include subhalo substructure
    double subhalo_m_floor = 1.0e7; // minimum clump mass
    double psi_min_fixed = -1.0; // fixed psi_min instead of m_floor/M, <= 0 off

    // substructure model. 1: reduced host, resolved clumps only above the
    //   dynamic r_thr floor. 3: model 1 plus the Gaussian unresolved term, exact
    //   in mean and variance below that floor. 4: brute, every subhalo resolved
    //   down to psi_min = m_floor/M with the host carved to M - sum_i m_i and no
    //   unresolved apparatus at all, costing ~45 ms/ray. 5: the same population
    //   as model 4, but rendering only the clumps that clear kappathr_sub, which
    //   is ~200x cheaper and sits at model 3 cost (DEFAULT).
    // Models 4 and 5 require subhalo_carve; model 1 requires it in this port.
    int subhalo_model = 5;
    // resolve down to m_floor with no dynamic floor. Meaningless for model 5,
    // which is brute by construction, and incompatible with model 3
    bool subhalo_brute = false;
    // mass conserving carve: build the host at M - sum_i m_i - M_u(r) from the
    // realized clump mass, so the total halo mass is M in every realization and
    // not merely in the mean. Intrinsic to models 4 and 5
    bool subhalo_carve = true;
    // refer psi to M_vir and extend the population to r_vir, per Jiang & van den
    // Bosch 2014. Gated to models 4 and 5 (see subhalo.h)
    bool subhalo_virial = true;
    // per-clump kappa threshold for model 5. subhalo_kappathr sets it absolutely
    // if positive, otherwise subhalo_kappathr_factor scales the host kappathr so
    // that it tracks the host rule. At 0.1 the clump count already falls ~1e4x
    // for a loss of ~0.1% in sigma, so the extra decade at 1.0 buys little
    double subhalo_kappathr = -1.0;
    double subhalo_kappathr_factor = 0.1;
    // model 3 only, ignored by models 4 and 5. Rescales the clump resolution
    // threshold. Clumps below it are folded into the analytic unresolved term,
    // which is exact in mean and variance but Gaussian, so the factor only
    // controls how much of the clump third cumulant is Gaussianized:
    // 0.16% at 1.0e-3, 0.73% at 1.0e-2
    double subhalo_factor = 1.0e-2;

    // ---- Anchor that enforces <kappa> = 0. The subtraction uses an empirical
    // mean, so one ray with kappa >> 1, already outside weak lensing validity,
    // shifts the whole batch by -2 kappa/n.
    // 0: running sum over encounters divided by Nreal, as the original code. It
    //   excludes the realized weak background, which is zero mean anyway.
    // 1: mean over the per-ray totals with kappa <= kappa_anchor_cut, so that no
    //   single ray can move the batch. The bias is of order cut/n (DEFAULT).
    // 2: use kappa_anchor_value directly, giving independent realizations.
    int kappa_anchor = 1;
    double kappa_anchor_cut = 1.0;
    double kappa_anchor_value = 0.0;

    // ---- Clustering bias of the discrete lens counts (ported from the emulator,
    // docs: bias_field_design_note / bias_window_design_plan / filament_bias_note).
    // bias_model 0 = legacy iid per-(z,M) log-normal cell modulation (DEFAULT,
    //   bit-identical to the original code); 1 = correlated 1D pencil-beam field
    //   delta_1D(chi) along the line of sight (segment-averaged, Cholesky-realized).
    // Set to 0 to recover the legacy iid layer bit-for-bit, and then also clear
    // bias_window, bias_weak and fil_bias below, which need bias_model = 1.
    int bias_model = 1;
    // Comoving transverse smoothing radius R_s of the field, kpc (bias_model = 1
    // only). The earlier 8441 = R_L(1e14) is a Lagrangian radius, not a
    // smoothing scale.
    double bias_Rperp = 20000.0;
    // Field smoothing window: 0 = transverse disk on k_perp (legacy prototype),
    // 1 = spherical top-hat, 2 = Gaussian (1/2 act on |k|). Nonzero requires
    // bias_model = 1 (throws).
    int bias_window = 1;
    // Weak (sub-threshold) arm drawn CONDITIONALLY on the same realized field
    // (Cox split of Campbell's theorem) instead of the unconditional Gaussian.
    // Requires bias_model = 1 (throws otherwise); no-op when bias = 0. Only NFW
    // field halos feed the weak arm; filaments carry none. At R_s = 20 Mpc it
    // carries much of the clustering effect, so a counts-only reading of the
    // field alone understates it.
    bool bias_weak = true;
    // Filaments ride the filament bias b_fil = filbias (PBS of pFCfil, q = 0.7)
    // instead of the halo bias. Requires bias_model = 1 (no-op in the iid layer);
    // reuses the same realized field, no new RNG draw.
    bool fil_bias = true;

    // probability distribution of lnmu, {lnmu, dP/dlnmu}
    vector<vector<double> > Plnmuf(cosmology &C, double zs, rgen &mt, int fil, int bias, int ell, int write);
    
    // MCMC likelihood analysis of the Hubble diagram
    void Hubble_diagram_fit(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int Ns, int Nburnin, int lens, int dm, rgen &mt, fs::path filename);
    
private:
    Subhalo S;
    
    // loglikelihood of the Hubble digram data
    double loglikelihood(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &par, int lens, int dm, rgen &mt);
    
};
