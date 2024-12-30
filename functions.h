#include <complex>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;

typedef mt19937_64 rgen;

const double PI = 3.141592653589793238463;

string to_string_prec(const double a, const int n);
double randomreal(double x1, double x2, rgen &mt);
double interpolate(double x, vector<vector<double> > &y);
double findrootG(double y, double dx, vector<vector<double> > &list);

vector<vector<double> > comovingdistance(function<double(double)> Hz, const int Nz, const double zmin, const double zmax);

class cosmology {
    public:
    
    // cosmological parameters
    double OmegaM;
    double OmegaB;
    double zeq;
    double sigma8;
    double h;
    double T0;
    double ns;
    double OmegaR;
    double OmegaL;
    double OmegaC;
    double H0;
    double rhoc;
    double rhoM0;
    
    // z dependencies
    double Az(double z) {
        return OmegaM*pow(1+z,3.0) + OmegaR*pow(1+z,4.0) + OmegaL;
    }
    double Hz(double z) {
        return H0*sqrt(Az(z));
    }
    double OmegaMz(double z) {
        return OmegaM*pow(1+z,3.0)/Az(z);
    }
    double OmegaRz(double z) {
        return OmegaR*pow(1+z,4.0)/Az(z);
    }
    double OmegaLz(double z) {
        return OmegaL/Az(z);
    }
    
    // growth function
    double Dg(double z) {
        return 5.0/2.0*OmegaMz(z)/(pow(OmegaMz(z),4.0/7.0) - OmegaLz(z) + (1+OmegaMz(z)/2.0)*(1+OmegaLz(z)/70.0))/(1+z)/0.7869370293916;
    }
    double delta(double z) {
        return 3.0*pow(12.0*PI,2.0/3.0)/20.0*(1+0.123*log10(OmegaL*pow(1+z,3.0)/(OmegaL*pow(1+z,3.0)+1-OmegaL)));
    }
    double deltac(double z) {
        return 3.0/5.0*pow(3.0*PI/2.0,2.0/3.0)/Dg(z);
    }
    
    // matter transfer function (astro-ph/9709112)
    double TM (double k);
    
    // matter power spectrum
    double Pk(double k, double deltaH) {
        return pow(306.535*k/H0,3.0+ns)*pow(deltaH*TM(k),2.0)*2.0*pow(PI,2.0)/pow(k,3.0);
    }
    double Deltak(double k, double deltaH) {
        return sqrt(pow(306.535*k/H0,3.0+ns)*pow(deltaH*TM(k),2.0));
    }
    
    // window function
    double W(double x) {
        return 3.0*(x*cos(x)-sin(x))/pow(x,3.0);
    }
    
    // variance of the matter fluctuations
    double sigma(double RM, double deltaH, int Nk);
    
    // variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    vector<vector<double> > sigmalist(int Nk, int NM, double Mmin, double Mmax);
    
    // Seth-Tormen HMF, {z,M,dndlnM}
    vector<vector<vector<double> > > hmflist(vector<vector<double> > &sigma3, int Nz, double zmin, double zmax);
};

