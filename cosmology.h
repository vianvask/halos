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
    double fB;
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
    
    // spherical collapse threshold
    double deltac(double z) {
        return 3.0/5.0*pow(3.0*PI/2.0,2.0/3.0)/Dg(z);
    }
    
    // elliptical collapse threshold (Eq. (1) of MNRAS 329, 61–75 (2002)), S>0
    double deltaell(double z, double S) {
        double alpha = 0.615;
        double beta = 0.485;
        double a = 0.707;
        double delta = deltac(z);

        return sqrt(a)*delta*(1.0 + beta*pow(a*pow(delta,2.0)/S,-alpha));
    }
    
private:
    
    // CDM matter transfer function (astro-ph/9709112)
    double TM(double k);
    
    // FDM matter transfer function (astro-ph/0003365)
    double kJ(double m22) {
        return 69.1e-3*pow(OmegaM*h*h/0.14/(1+zeq),1.0/4.0)*sqrt(m22);
    }
    double km22(double m22) {
        return 1.0/1.61*pow(m22,-1.0/18.0)*kJ(m22);
    }
    double TF(double k, double m22) {
        return cos(pow(k/km22(m22),3.0))/(1+pow(k/km22(m22),8.0));
    }
    double TMF(double k, double m22) {
        return TF(k,m22)*TM(k);
    }
    
    // WDM matter transfer function (astro-ph/0501562)
    double km3(double m3) {
        return 1.0/(49.0/h*pow(OmegaC/0.25,0.11)*pow(h/0.7,1.22)*pow(m3,-1.11));
    }
    double TW(double k, double m3) {
        double mu = 1.12;
        return pow(1.0+pow(k/km3(m3),2.0*mu),-5.0/mu);
        
    }
    double TMW(double k, double m3) {
        return TW(k,m3)*TM(k);
        
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
    
    // real space top hat window function
    double W(double x) {
        return 3.0*(x*cos(x)-sin(x))/pow(x,3.0);
    }
    double DW(double x) {
        return -9.0*cos(x)/pow(x,3.0) + 9.0*sin(x)/pow(x,4.0) - 3.0*sin(x)/pow(x,2.0);
    }
    
    // smooth k space window function
    double Ws(double x, double c, double b) {
        return 1.0/(1.0+pow(x/c,b));
    }
    double DWs(double x, double c, double b) {
        return -b*pow(x/c,b)/(x*pow(1.0+pow(x/c,b),2.0));
    }
    
    // variance of the matter fluctuations, m22: FDM mass in 10^-22 eV, m3: WDM mass in keV
    vector<double> sigmaC(double M, double deltaH);
    vector<double> sigmaF(double M, double deltaH, double m22);
    vector<double> sigmaW(double M, double deltaH, double m3);

    // variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    vector<vector<double> > sigmalistf(double m22, double m3);
    
    vector<double> zlistf();
    vector<double> m22listf();
    vector<double> m3listf();
    
    // CDM halo consentration parameter (1601.02624)
    double cons(double M, double z);
    vector<vector<vector<double> > > conslistf();
    vector<vector<vector<double> > > conslist;
    
    // NFW scale radius and density and their derivatives, {z, M, r_s, dr_s/dM, rho_s, drho_s/dM}
    vector<vector<vector<double> > > NFWlistf();
    
    // first crossing probability with ellipsoidal collapse
    double pFC(double delta, double S);
    double pFC(double deltap, double Sp, double delta, double S);

    // halo mass function, {z,M,dndlnM}
    vector<vector<vector<double> > > HMFlistf();

    // growth rate of the halo through mergers with smaller halos, {z,M,dM/dt}
    vector<vector<vector<double> > > dotMlistf();
    
    vector<double> growbymergers(vector<double> &MJ, double z, double zp);

    vector<vector<double> > dclist();
    vector<vector<double> > tlist();
    
    vector<vector<double> > zdc;
    vector<vector<double> > zt;
    
public:
    // star formation rate
    double fstar(double z, double M, double Mc, double Mt, double epsilon, double alpha, double beta);

    // derivative of the star formation rate, df_*/dM
    double Dfstarperfstar(double M, double Mc, double Mt, double epsilon, double alpha, double beta);
    
    vector<vector<double> > evolvestellarmass(double Mc, double epsilon, double alpha, double beta);
    vector<vector<double> > evolveBHmass(double Mc, double epsilon, double alpha, double beta, double fEdd, double facc1, double facc2);
    
    // list of redshifts, same as for HMF and dotM lists
    vector<double> zlist;
    vector<vector<double> > sigmalist;
    vector<vector<vector<double> > > HMFlist;
    vector<vector<vector<double> > > dotMlist;
    vector<vector<vector<double> > > NFWlist;
    
    // list of redshifts, same as for P(lnmu,z) list
    vector<double> Zlist;
    vector<vector<vector<double> > > Plnmuz;
    
    // list of FDM masses, same as for HMF_FDM and dotM_FDM lists
    vector<double> m22list;
    vector<vector<vector<double> > > FDMsigmalist;
    vector<vector<vector<vector<double> > > > FDMHMFlist;
    vector<vector<vector<vector<double> > > > FDMdotMlist;
    
    // list of WDM masses, same as for HMF_WDM and dotM_WDM lists
    vector<double> m3list;
    vector<vector<vector<double> > > WDMsigmalist;
    vector<vector<vector<vector<double> > > > WDMHMFlist;
    vector<vector<vector<vector<double> > > > WDMdotMlist;
    
    void initialize() {
        OmegaR = OmegaM/(1+zeq);
        OmegaL = 1.0 - OmegaM - OmegaR;
        OmegaC = OmegaM - OmegaB;
        fB = OmegaB/OmegaM;
        
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
        deltaH8 = sigma8/sigmaC(M8, 1.0)[0];
        
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
            deltaH8 = sigma8/sigmaF(M8, 1.0, m22)[0];
                        
            // halo mass function and halo growth rate
            sigmalist = sigmalistf(m22, 0.0);
            HMFlist = HMFlistf();
            dotMlist = dotMlistf();
            
            FDMsigmalist.push_back(sigmalist);
            FDMHMFlist.push_back(HMFlist);
            FDMdotMlist.push_back(dotMlist);
        }
    }
    
    void WDM_halos() {
        for (double m3 : m3list) {
            // fix deltaH to match the input sigma8
            deltaH8 = sigma8/sigmaW(M8, 1.0, m3)[0];
                        
            // halo mass function and halo growth rate
            sigmalist = sigmalistf(0.0, m3);
            HMFlist = HMFlistf();
            dotMlist = dotMlistf();
            
            WDMsigmalist.push_back(sigmalist);
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
