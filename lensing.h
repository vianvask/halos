#include "cosmology.h"

// lensing amplification lnmu and its derivative dlnmu/dlnM
vector<double> kappa(cosmology &C, double zs, double zl, double r, double M, double rS);

// probability distribution {lnmu, P1(lnmu)}
vector<vector<double> > dNdkappa(cosmology &C, int N, double zs, double lnmumax, double rS);
