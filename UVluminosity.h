
// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
vector<vector<vector<double> > > PhiUV(cosmology &C, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma, double zc);
vector<vector<vector<double> > > PhiUV(cosmology &C, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma, double zbreak, double m);

vector<double> UFfit(cosmology &C, vector<vector<double> > &priors, vector<double> &steps, int Nsteps, int Nbi, int Nchains);
