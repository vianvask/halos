
// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
vector<vector<vector<double> > > PhiUV(cosmology &C, double logMt, double Mc, double epsilon, double alpha, double beta, double gamma, double zc, double fkappa, double ze, double z0, double sigmaUV, double logm);

vector<double> UFfit(cosmology &C, vector<vector<double> > &priors, vector<double> &steps, int Nsteps, int Nbi, int Nchains);

double loglikelihood(cosmology &C, vector<double> &params);
