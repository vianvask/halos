#include "basics.h"

vector<vector<double> > comovingdistance(function<double(double)> Hz, const int Nz, const double zmin, const double zmax);

class cosmology {

public:
    int Nk, NM, Nz;
    double Mmin, Mmax, zmin, zmax;
    
    int Nm22 = 0;
    double m22min, m22max;
    
    int Nm3 = 0;
    double m3min, m3max;
    
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
    
    double M8;
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
    
    // k scale for FDM
    double km22(double m22) {
        return 1.0/1.61*pow(m22,-1.0/18.0)*kJ(m22);
    }
    
    // k scale for WDM
    double km3(double m3) {
        return 1.0/72.7*pow(m3,1.11);
    }
    
    // FDM matter transfer function (astro-ph/0003365)
    double TMF (double k, double m22) {
        return cos(pow(k/km22(m22),3.0))/(1+pow(k/km22(m22),8.0))*TM(k);
    }
    
    // WDM matter transfer function
    double TMW (double k, double m3) {
        double mu = 1.12;
        return pow(1.0+pow(k/km3(m3),2.0*mu),-5.0/mu)*TM(k);
        
    }
    
    // CDM matter power spectrum
    double Deltak(double k, double deltaH) {
        return sqrt(pow(306.535*k/H0,3.0+ns)*pow(deltaH*TM(k),2.0));
    }
    
    // FDM matter power spectrum
    double DeltakF(double k, double deltaH, double m22) {
        return sqrt(pow(306.535*k/H0,3.0+ns)*pow(deltaH*TMF(k,m22),2.0));
    }
    
    // WDM matter power spectrum
    double DeltakW(double k, double deltaH, double m3) {
        return sqrt(pow(306.535*k/H0,3.0+ns)*pow(deltaH*TMW(k,m3),2.0));
    }
    
    // CDM window function
    double W(double x) {
        return 3.0*(x*cos(x)-sin(x))/pow(x,3.0);
    }
    
    // FDM window function
    double WF(double x) {
        double a = 0.43;
        double b1 = 7.0;
        double b2 = 1.0;
        return 1.0-1.0/pow(1.0+pow(a*x,-b1/b2),b2);
    }
    
    // WDM window function
    double WW(double x) {
        double a = 0.43;
        double b1 = 7.0;
        double b2 = 1.0;
        return 1.0-1.0/pow(1.0+pow(a*x,-b1/b2),b2);
    }
    
    // variance of the matter fluctuations, m22: FDM mass in 10^-22 eV, m3: WDM mass in keV
    double sigmaf(double M, double deltaH, double m22, double m3);
    
    // variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    vector<vector<double> > sigmalistf(double m22, double m3);
    vector<vector<double> > sigmalist;
    
    vector<double> zlistf();
    vector<double> m22listf();
    vector<double> m3listf();
    
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
    
    // list of redshifts, same as for HMF and dotM lists
    vector<double> zlist;
    vector<vector<vector<double> > > HMFlist;
    vector<vector<vector<double> > > dotMlist;
    vector<vector<vector<double> > > NFWlist;
    
    // list of redshifts, same as for P(lnmu,z) list
    vector<double> Zlist;
    vector<vector<vector<double> > > Plnmuz;
    
    // list of FDM masses, same as for HMF_FDM and dotM_FDM lists
    vector<double> m22list;
    vector<vector<vector<vector<double> > > > FDMHMFlist;
    vector<vector<vector<vector<double> > > > FDMdotMlist;
    
    // list of WDM masses, same as for HMF_WDM and dotM_WDM lists
    vector<double> m3list;
    vector<vector<vector<vector<double> > > > WDMHMFlist;
    vector<vector<vector<vector<double> > > > WDMdotMlist;
    
    void initialize() {
        OmegaR = OmegaM/(1+zeq);
        OmegaL = 1.0 - OmegaM - OmegaR;
        OmegaC = OmegaM - OmegaB;
        H0 = 0.000102247*h;
        rhoc = 277.394*pow(h,2.0);
        rhoM0 = OmegaM*rhoc;
        M8 = 4.0*PI/3.0*pow(8000.0/h,3.0)*rhoM0;
        
        zlist = zlistf();
        m22list = m22listf();
        m3list = m3listf();
        
        zdc = dclist();
        zt = tlist();
    }
       
    void CDM_halos() {
        // fix deltaH to match the input sigma8
        deltaH8 = sigma8/sigmaf(M8, 1.0, 0.0, 0.0);
        
        // halo mass function and halo growth rate
        sigmalist = sigmalistf(0.0, 0.0);
        HMFlist = HMFlistf();
        dotMlist = dotMlistf();

        // NFW halo parameters
        conslist = conslistf();
        NFWlist = NFWlistf();
    }
    
    void FDM_halos() {
        for (double m22 : m22list) {
            // fix deltaH to match the input sigma8
            deltaH8 = sigma8/sigmaf(M8, 1.0, m22, 0.0);
                        
            // halo mass function and halo growth rate
            sigmalist = sigmalistf(m22, 0.0);
            HMFlist = HMFlistf();
            dotMlist = dotMlistf();
            
            FDMHMFlist.push_back(HMFlist);
            FDMdotMlist.push_back(dotMlist);
        }
    }
    
    void WDM_halos() {
        for (double m3 : m3list) {
            // fix deltaH to match the input sigma8
            deltaH8 = sigma8/sigmaf(M8, 1.0, 0.0, m3);
                        
            // halo mass function and halo growth rate
            sigmalist = sigmalistf(0.0, m3);
            HMFlist = HMFlistf();
            dotMlist = dotMlistf();
            
            WDMHMFlist.push_back(HMFlist);
            WDMdotMlist.push_back(dotMlist);
        }
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
