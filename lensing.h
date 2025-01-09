#include "cosmology.h"

// lensing amplification {lnmu, dlnmu/dlnM}
vector<double> lnmu(cosmology C, double zs, double zl, double r, double M);

// probability distribution {lnmu, P1(lnmu)}
vector<vector<double> > P1(int N);

