#include "basics.h"

vector<vector<double> > comovingdistance(function<double(double)> Hz, const int Nz, const double zmin, const double zmax);

class cosmology {

public:
    int Nk = 1000;
    int NM, Nz;
    double Mmin, Mmax, zmin, zmax;
    
    int Nm22 = 0;
    double m22min, m22max;
    
    int Nm3 = 0;
    double m3min, m3max;
    
    int Nkc = 0;
    double kcmin, kcmax;
    
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
    
    // FDM matter transfer function (2201.10238, old from astro-ph/0003365)
    double kJ(double m22) {
        return 66.5e-3*pow(OmegaC*h*h/0.12/(1+zeq),1.0/4.0)*sqrt(m22);
    }
    double km22(double m22) {
        double A = 2.22*pow(m22,1.0/25.0 - 0.001*log(m22));
        return kJ(m22)/A;
        //return 1.0/1.61*pow(m22,-1.0/18.0)*kJ(m22);
    }
    double TF(double k, double m22) {
        double n = 5.0/2.0;
        double x = k/km22(m22);
        double B = 0.16*pow(m22,-1.0/20.0);
        return sin(pow(x,n))/pow(x,n)/(1+B*pow(x,6-n));
        //return cos(pow(k/km22(m22),3.0))/(1+pow(k/km22(m22),8.0));
    }
    double TMF(double k, double m22) {
        return TF(k,m22)*TM(k);
    }
    
    // WDM matter transfer function (astro-ph/0501562)
    double km3(double m3) {
        double fix = 3.3*0.43; // fixes the problem with M_hm (see 1801.02547)
        return 1.0/(fix*49.0/h*pow(OmegaC/0.25,0.11)*pow(h/0.7,1.22)*pow(m3,-1.11));
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
    
    // white noise enhanced matter power spectrum
    double DeltakE(double k, double deltaH, double kc) {
        return Deltak(k, deltaH) + pow(k/kc,3.0)*Deltak(kc, deltaH);
    }
    
    // real space top hat window function
    double W(double x) {
        return 3.0*(x*cos(x)-sin(x))/pow(x,3.0);
    }
    double DW(double x) {
        return -9.0*cos(x)/pow(x,3.0) + 9.0*sin(x)/pow(x,4.0) - 3.0*sin(x)/pow(x,2.0);
    }
    
    // smooth k space window function
    double Ws(double x) {
        double c = 1.0/0.43;
        double b = 6.0;
        return 1.0/(1.0+pow(x/c,b));
    }
    double DWs(double x) {
        double c = 1.0/0.43;
        double b = 6.0;
        return -b*pow(x/c,b)/(x*pow(1.0+pow(x/c,b),2.0));
    }
    
    // variance of the matter fluctuations, m22: FDM mass in 10^-22 eV, m3: WDM mass in keV
    vector<double> sigmaC(double M, double deltaH);
    vector<double> sigmaF(double M, double deltaH, double m22);
    vector<double> sigmaW(double M, double deltaH, double m3);
    vector<double> sigmaE(double M, double deltaH, double kc);

    // variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    vector<vector<double> > sigmalistf(double m22, double m3, double kc);
    
    // generates a list in log scale
    vector<double> loglist(double xmin, double xmax, int Nx);

    // CDM halo consentration parameter (1601.02624), {jz,jM} -> {c, dc/dM}
    double cons(double M, double z);
    vector<vector<vector<double> > > conslistf();
    vector<vector<vector<double> > > conslist;
    
    // NFW scale radius and density and their derivatives, {jz,jM} -> {r_s, dr_s/dM, rho_s, drho_s/dM}
    vector<vector<vector<double> > > NFWlistf();
    
    // first crossing probability with ellipsoidal collapse
    double pFC(double delta, double S);
    double pFC(double deltap, double Sp, double delta, double S);

    // halo mass function and growth rate, {jz,jM} -> {dn/dlnM, dotM, dotM/dM}
    vector<vector<vector<double> > > HMFlistf();
    
    vector<vector<double> > dclist();
    vector<vector<double> > tlist();
    
    vector<vector<double> > zdc;
    vector<vector<double> > zt;
    
public:
    // star formation rate
    double fstar(double M, double Mc, double Mt, double epsilon, double alpha, double beta);

    // derivative of the star formation rate, df_*/dM
    double Dfstarperfstar(double M, double Mc, double Mt, double epsilon, double alpha, double beta);
    
    vector<vector<double> > evolvestellarmass(double Mc, double Mt, double epsilon, double alpha, double beta);
    vector<vector<double> > evolveBHmass(double Mc, double Mt, double epsilon, double alpha, double beta, double fEdd, double facc1, double facc2);
    
    // list of redshifts and halo masses
    vector<double> zlist;
    vector<double>  Mlist;
    
    vector<vector<double> > sigmalist;
    vector<vector<vector<double> > > HMFlist;
    vector<vector<vector<double> > > NFWlist;
    
    // list of redshifts, same as for P(lnmu,z) list
    vector<double> Zlist;
    vector<vector<vector<double> > > Plnmuz;
    
    vector<double> m22list;
    vector<vector<vector<double> > > FDMsigmalist;
    vector<vector<vector<vector<double> > > > FDMHMFlist;
    
    vector<double> m3list;
    vector<vector<vector<double> > > WDMsigmalist;
    vector<vector<vector<vector<double> > > > WDMHMFlist;
    
    vector<double> kclist;
    vector<vector<vector<double> > > EDMsigmalist;
    vector<vector<vector<vector<double> > > > EDMHMFlist;
       
    void CDM_halos() {
        // fix deltaH to match the input sigma8
        deltaH8 = sigma8/sigmaC(M8, 1.0)[0];
        
        // halo mass function and halo growth rate
        sigmalist = sigmalistf(0.0, 0.0, 0.0);
        HMFlist = HMFlistf();
        
        writeToFile(sigmalist, "sigma_CDM.dat");
        writeToFile(zlist, Mlist, HMFlist, "HMF_CDM.dat");
        
        // NFW halo parameters
        conslist = conslistf();
        NFWlist = NFWlistf();
    }
    
    void FDM_halos() {
        for (double m22 : m22list) {
            // fix deltaH to match the input sigma8
            deltaH8 = sigma8/sigmaF(M8, 1.0, m22)[0];
                        
            // halo mass function and halo growth rate
            sigmalist = sigmalistf(m22, 0.0, 0.0);
            HMFlist = HMFlistf();
            
            FDMsigmalist.push_back(sigmalist);
            FDMHMFlist.push_back(HMFlist);
        }
        writeToFile(m22list, FDMsigmalist, "sigma_FDM.dat");
        writeToFile(m22list, zlist, Mlist, FDMHMFlist, "HMF_FDM.dat");
    }
    
    void WDM_halos() {
        for (double m3 : m3list) {
            // fix deltaH to match the input sigma8
            deltaH8 = sigma8/sigmaW(M8, 1.0, m3)[0];
                        
            // halo mass function and halo growth rate
            sigmalist = sigmalistf(0.0, m3, 0.0);
            HMFlist = HMFlistf();
            
            WDMsigmalist.push_back(sigmalist);
            WDMHMFlist.push_back(HMFlist);
        }
        
        writeToFile(m3list, WDMsigmalist, "sigma_WDM.dat");
        writeToFile(m3list, zlist, Mlist, WDMHMFlist, "HMF_WDM.dat");
    }
    
    void EDM_halos() {
        for (double kc : kclist) {
            // fix deltaH to match the input sigma8
            deltaH8 = sigma8/sigmaE(M8, 1.0, kc)[0];
                        
            // halo mass function and halo growth rate
            sigmalist = sigmalistf(0.0, 0.0, kc);
            HMFlist = HMFlistf();
            
            EDMsigmalist.push_back(sigmalist);
            EDMHMFlist.push_back(HMFlist);
        }
        
        writeToFile(kclist, EDMsigmalist, "sigma_EDM.dat");
        writeToFile(kclist, zlist, Mlist, EDMHMFlist, "HMF_EDM.dat");
    }
    
    void initialize(int dm) {
        OmegaR = OmegaM/(1+zeq);
        OmegaL = 1.0 - OmegaM - OmegaR;
        OmegaC = OmegaM - OmegaB;
        fB = OmegaB/OmegaM;
        
        H0 = 0.000102247*h;
        rhoc = 277.394*pow(h,2.0);
        rhoM0 = OmegaM*rhoc;
        M8 = 4.0*PI/3.0*pow(8000.0/h,3.0)*rhoM0;
        
        zlist = loglist(zmin,zmax,Nz);
        Mlist = loglist(Mmin,Mmax,NM);
                
        zdc = dclist();
        zt = tlist();
        
        if (dm == 0) {
            CDM_halos();
        }
        if (dm == 1) {
            m22list = loglist(m22min,m22max,Nm22);
            FDM_halos();
        }
        if (dm == 2) {
            m3list = loglist(m3min,m3max,Nm3);
            WDM_halos();
        }
        if (dm == 3) {
            kclist = loglist(kcmin,kcmax,Nkc);
            EDM_halos();
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
