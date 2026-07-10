class cosmology;

class Subhalo {
    
public:
    double alpha = -0.82, beta = 50.0, omega = 4.0, psi_res = 1.0e-4, psi_max = 0.1;
    double m_floor = 1.0e7;
    int Nu = 128, NyW = 48;
    
    vector<vector<double> > gnorm;
    vector<vector<double> > r200h;
    vector<vector<double> > chost;
    vector<vector<double> > r_thr;
    vector<vector<vector<double> > > invRad;
    vector<vector<vector<double> > > muW;
    vector<vector<vector<double> > > sW;
    vector<vector<array<double,2> > > lyW;
    vector<vector<double> > fsb;
    double log_Mmin = 0.0;
    double inv_dlogM = 0.0;
    
    void precompute(cosmology &C, double zs, double kappathr, double kappathr_host);
    void buildWsubBin(cosmology &C, double zs, int jz, int jM, double kappathr_host);
    void wsubTerm(int jz, int jM, double r, double &mu, double &sigma);
    int addClumps(cosmology &C, int jz, int jM, double zl, double M, double Sigmac, double r, double phi, rgen &mt, double &kappa, double &gamma1, double &gamma2);
    
};
