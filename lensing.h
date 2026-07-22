
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
    
    // probability distribution of lnmu, {lnmu, dP/dlnmu}
    vector<vector<double> > Plnmuf(cosmology &C, double zs, rgen &mt, int fil, int bias, int ell, int write);
    
    // MCMC likelihood analysis of the Hubble diagram
    void Hubble_diagram_fit(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int Ns, int Nburnin, int lens, int dm, rgen &mt, fs::path filename);
    
private:
    Subhalo S;
    
    // loglikelihood of the Hubble digram data
    double loglikelihood(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &par, int lens, int dm, rgen &mt);
    
};
