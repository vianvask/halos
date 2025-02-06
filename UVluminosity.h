
// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
vector<vector<vector<double> > > PhiUV(cosmology &C, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma, double zbreak);
vector<vector<vector<double> > > PhiUV(cosmology &C, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma, double zbreak, double m22);

// read UV luminosity data, {MUV, Phi, +sigmaPhi, -sigmaPhi}
vector<vector<double> > readUVdata(vector<string> filenames);

// Metropolis-Hastings MCMC sampler of the UV luminosity fit likelihood
vector<vector<double> > mcmc_sampling(cosmology &C, vector<vector<double> > &data, vector<double> &initial, vector<double> &steps , vector<vector<double> > &priors, int Ns, rgen &mt);

double loglikelihood(cosmology &C, vector<vector<double> > &data, vector<double> &params);
