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
    double sigmaf(double M, double deltaH);
    
    // variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    vector<vector<double> > sigmalistf();
    vector<vector<double> > sigmalist;
    
    // halo consentration parameter (1601.02624)
    double cons(double M, double z);
    vector<vector<vector<double> > > conslistf();
    vector<vector<vector<double> > > conslist;
    
    // NFW scale radius and density and their derivatives, {z, M, r_s, dr_s/dM, rho_s, drho_s/dM}
    vector<vector<vector<double> > > NFWlistf();
    
    // Seth-Tormen HMF, {z,M,dndlnM}
    vector<vector<vector<double> > > hmflistf();
    
    double dPpdM(vector<double> &sigma, double z, vector<double> &sigmap, double zp);
    double DeltaM(vector<double> &sigmap, double z, double dz);
    
    // growth rate of the halo through mergers with smaller halos, {z,M,dM/dt}
    vector<vector<vector<double> > > dotMlistf();
    
    // star formation rate
    double fstar(double M, double Mc, double epsilon, double alpha, double beta) {
        return epsilon/(pow(M/Mc,alpha) + pow(M/Mc,beta));
    }
    // derivative of the star formation rate, df_*/dM
    double Dfstar(double M, double Mc, double epsilon, double alpha, double beta) {
        return -(alpha*pow(M/Mc,alpha) + beta*pow(M/Mc,beta))*pow(fstar(M,Mc,epsilon,alpha,beta),2.0)/(epsilon*M);
    }
    
    double kappaUV = 1.15e-22; // Msun s erg^-1 Myr^-1
    
    // the UV magnitude
    double MUV(double M, double dotM, double Mc, double epsilon, double alpha, double beta) {
        return 51.63 - 1.08574*log(fstar(M,Mc,epsilon,alpha,beta)*dotM/kappaUV);
    }
        
    // derivative of the UV magnitude, dM_UV/dM
    double DMUV(double M, double dotM, double DdotM, double Mc, double epsilon, double alpha, double beta) {
        return -1.08574*(DdotM/dotM + Dfstar(M,Mc,epsilon,alpha,beta)/fstar(M,Mc,epsilon,alpha,beta));
    }
    
    // dust extinction, MUV -> MUV - AUV (see 1406.1503 and Table 3 of 1306.2950)
    double AUV(double MUV, double z) {
        double C0 = 4.4, C1 = 2.0, sigmabeta = 0.34;
        double beta = -exp(0.17*(19.5+MUV)/(1.54+0.075*z))*(1.54+0.075*z);
        return max(0.0, C0 + 0.2*log(10)*pow(C1*sigmabeta,2.0) + C1*beta);
    }

    // UV luminosity function as a function of halo mass M (see e.g. 1906.06296)
    double UVLF(double M, double dotM, double DdotM, double dndlnM, double Mt, double Mc, double epsilon, double alpha, double beta) {
        return -exp(-Mt/M)*dndlnM/M/DMUV(M,dotM,DdotM,Mc,epsilon,alpha,beta);
    }
    
    vector<vector<double> > dclist();
    vector<vector<double> > tlist();
    
    vector<vector<double> > zdc;
    vector<vector<double> > zt;
    
public:
    
    vector<vector<vector<double> > > hmflist;
    vector<vector<vector<double> > > dotMlist;
    vector<vector<vector<double> > > NFWlist;
    
    // UV luminosity function, {z,M,MUV,AUV,Phi}
    vector<vector<vector<double> > > UVLFlistf(double Mt, double Mc, double epsilon, double alpha, double beta);
    
    void initialize() {
        OmegaR = OmegaM/(1+zeq);
        OmegaL = 1.0 - OmegaM - OmegaR;
        OmegaC = OmegaM - OmegaB;
        H0 = 0.000102247*h;
        rhoc = 277.394*pow(h,2.0);
        rhoM0 = OmegaM*rhoc;
        
        // fix deltaH to match the input sigma8
        double M8 = 4.0*PI/3.0*pow(8000.0/h,3.0)*rhoM0;
        deltaH8 = sigma8/sigmaf(M8, 1.0);
        
        zdc = dclist();
        zt = tlist();
                
        sigmalist = sigmalistf();
        conslist = conslistf();
        NFWlist = NFWlistf();
        hmflist = hmflistf();
        dotMlist = dotMlistf();
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
