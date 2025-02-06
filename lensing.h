
// NFW lensing amplification and its derivative, {kappa, dkappa/dr}
vector<double> kappa(cosmology &C, double zs, double zl, double r, double M, double rS);

// single lens probability distribution normalized to number of haloes, {kappa, dN/dkappa}
vector<vector<double> > dNdkappa(cosmology &C, int N, double zs, double lnmumax, double rS);

// probability distribution of lnmu, {lnmu, dP/dlnmu}
vector<vector<double> > Plnmuf(cosmology &C, int Nx, double zs, double kappamax, double rS, int Nreal, int Nbins, rgen &mt);

// read or generate lensing amplification distribution
vector<vector<vector<double> > > getPlnmu(cosmology &C, double rS, int Nkappa, int Nreal, int Nbins);

