// probability distribution of lnmu, {lnmu, dP/dlnmu}
vector<vector<double> > Plnmuf(cosmology &C, double zs, int Nhalos, int Nreal, int Nbins, rgen &mt, int write);

// MCMC likelihood analysis of the Hubble diagram
void Hubble_diagram_fit(cosmology &C, double DLthr, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int Ns, int Nburnin, int lensing, int dm, rgen &mt, string filename);
