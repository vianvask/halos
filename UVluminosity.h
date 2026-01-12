
vector<vector<vector<double> > > getPlnmu(fs::path filename);

// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
void writeUVLF(cosmology &C, double logMt, double Mc, double epsilon, double alpha, double beta, double gamma, double zc, double fkappa, double ze, double z0, double sigmaUV, double logm, int dm);

vector<double> UVLFfit(cosmology &C, vector<vector<double> > &priors, int Nsteps, int Nburnin, int Nchains, double xstep, int dm);
