#include "cosmology.h"

// lensing amplification lnmu and its derivative dlnmu/dlnM
vector<double> lnmu(cosmology &C, double zs, double zl, double r, double M);

// probability distribution {lnmu, P1(lnmu)}
vector<vector<double> > dNdlnmu(cosmology &C, int N, double zs, double lnmuthr, double lnmumax);
