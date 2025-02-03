#include "basics.h"

vector<vector<double> > comovingdistance(function<double(double)> Hz, const int Nz, const double zmin, const double zmax);

class cosmology {

public:
    int Nk, NM, Nz;
    double Mmin, Mmax, zmin, zmax;
    
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
    
    double deltaH8;
    
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
    
private:
    
    // CDM matter transfer function (astro-ph/9709112)
    double TM (double k);
    
    // Jeans length
    double kJ(double m22) {
        return 66.5e-3/pow(1+zeq,1.0/4.0)*sqrt(m22);
    }
    
    // FDM matter transfer function (astro-ph/0003365)
    double TMF (double k, double m22);
    
    // CDM matter power spectrum
    double Deltak(double k, double deltaH) {
        return sqrt(pow(306.535*k/H0,3.0+ns)*pow(deltaH*TM(k),2.0));
    }
    
    // FDM matter power spectrum
    double Deltak(double k, double deltaH, double m22) {
        return sqrt(pow(306.535*k/H0,3.0+ns)*pow(deltaH*TMF(k,m22),2.0));
    }
    
    // CDM window function
    double W(double x) {
        return 3.0*(x*cos(x)-sin(x))/pow(x,3.0);
    }
    
    // FDM window function
    double WF(double x) {
        double a = 0.447;
        double b = 5.0;
        return 1.0/(1.0+pow(a*x,b));
    }
    
    // variance of the matter fluctuations
    double sigmaf(double M, double deltaH, double m22);
    
    // variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    vector<vector<double> > sigmalistf(double m22);
    vector<vector<double> > sigmalist;
    
    vector<double> zlistf();
    
    // CDM halo consentration parameter (1601.02624)
    double cons(double M, double z);
    vector<vector<vector<double> > > conslistf();
    vector<vector<vector<double> > > conslist;
    
    // NFW scale radius and density and their derivatives, {z, M, r_s, dr_s/dM, rho_s, drho_s/dM}
    vector<vector<vector<double> > > NFWlistf();
    
    // Seth-Tormen HMF, {z,M,dndlnM}
    vector<vector<vector<double> > > HMFlistf();
    
    double dPpdM(vector<double> &sigma, double z, vector<double> &sigmap, double zp);
    double DeltaM(vector<double> &sigmap, double z, double dz);
    
    // growth rate of the halo through mergers with smaller halos, {z,M,dM/dt}
    vector<vector<vector<double> > > dotMlistf();
    
    vector<vector<double> > dclist();
    vector<vector<double> > tlist();
    
    vector<vector<double> > zdc;
    vector<vector<double> > zt;
    
public:
    
    vector<double> zlist;
    vector<vector<vector<double> > > HMFlist;
    vector<vector<vector<double> > > dotMlist;
    vector<vector<vector<double> > > NFWlist;
    
    void initialize(double m22) {
        OmegaR = OmegaM/(1+zeq);
        OmegaL = 1.0 - OmegaM - OmegaR;
        OmegaC = OmegaM - OmegaB;
        H0 = 0.000102247*h;
        rhoc = 277.394*pow(h,2.0);
        rhoM0 = OmegaM*rhoc;
        
        // fix deltaH to match the input sigma8
        double M8 = 4.0*PI/3.0*pow(8000.0/h,3.0)*rhoM0;
        deltaH8 = sigma8/sigmaf(M8, 1.0, m22);
                        
        zlist = zlistf();
        sigmalist = sigmalistf(m22);
        HMFlist = HMFlistf();
        dotMlist = dotMlistf();
        
        zdc = dclist();
        zt = tlist();
        conslist = conslistf();
        NFWlist = NFWlistf();
    }
    
    // luminosity distance
    double DL(double z) {
        return (1+z)*interpolate(z, zdc);
    }
    
    // age of the universe
    double age(double z) {
        return interpolate(z, zt);
    }
};
