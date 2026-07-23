
#include "subhalo.h"

class lensing {
    
public:
    int Nreal; // realizations
    int Nhalos; // number of halos in each realization
    int Nbins; // P(lnmu) bins
    bool subhalo = false; // include subhalo substructure
    double subhalo_m_floor = 1.0e7; // minimum clump mass
    // rescales the clump resolution threshold. Clumps below it are folded into the
    // analytic unresolved term, which is exact in mean and variance but Gaussian, so
    // the factor only controls how much of the clump third cumulant is Gaussianized:
    // 0.16% at 1.0e-3, 0.73% at 1.0e-2. Either value is adequate; 1.0e-3 is the default
    // because the threshold it scales grows with zs, so a fixed factor coarsens at high zs.
    double subhalo_factor = 1.0e-3;

    // ---- Clustering bias of the discrete lens counts (ported from the emulator,
    // docs: bias_field_design_note / bias_window_design_plan / filament_bias_note).
    // bias_model 0 = legacy iid per-(z,M) log-normal cell modulation (DEFAULT,
    //   bit-identical to the original code); 1 = correlated 1D pencil-beam field
    //   delta_1D(chi) along the line of sight (segment-averaged, Cholesky-realized).
    int bias_model = 0;
    // Comoving transverse smoothing radius R_s of the field, kpc (bias_model = 1
    // only). Production = 20000 (20 Mpc spherical top-hat).
    double bias_Rperp = 8441.0;
    // Field smoothing window: 0 = transverse disk on k_perp (legacy prototype),
    // 1 = spherical top-hat, 2 = Gaussian (1/2 act on |k|). Nonzero requires
    // bias_model = 1 (throws). Production = 1.
    int bias_window = 0;
    // Weak (sub-threshold) arm drawn CONDITIONALLY on the same realized field
    // (Cox split of Campbell's theorem) instead of the unconditional Gaussian.
    // Requires bias_model = 1 (throws otherwise); no-op when bias = 0. Production
    // = true. Only NFW field halos feed the weak arm; filaments carry none.
    bool bias_weak = false;
    // Filaments ride the filament bias b_fil = filbias (PBS of pFCfil, q = 0.7)
    // instead of the halo bias. Requires bias_model = 1 (no-op in the iid layer);
    // reuses the same realized field, no new RNG draw. Production = true.
    bool fil_bias = false;

    // probability distribution of lnmu, {lnmu, dP/dlnmu}
    vector<vector<double> > Plnmuf(cosmology &C, double zs, rgen &mt, int fil, int bias, int ell, int write);
    
    // MCMC likelihood analysis of the Hubble diagram
    void Hubble_diagram_fit(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int Ns, int Nburnin, int lens, int dm, rgen &mt, fs::path filename);
    
private:
    Subhalo S;
    
    // loglikelihood of the Hubble digram data
    double loglikelihood(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &par, int lens, int dm, rgen &mt);
    
};
