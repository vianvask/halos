class cosmology;

class Subhalo {
    
public:
    double alpha = -0.82, beta = 50.0, omega = 4.0, psi_res = 1.0e-4, psi_max = 0.1;
    double m_floor = 1.0e7;
    int Nu = 128;
    
    vector<vector<double> > gnorm;
    vector<vector<double> > Nsub;
    vector<vector<double> > r200h;
    vector<vector<double> > r_thr;
    vector<vector<vector<double> > > invRad;
    double log_Mmin = 0.0;
    double inv_dlogM = 0.0;
    
    void precompute(cosmology &C, double zs, double kappathr);
    double resolvedFraction(cosmology &C, int jz, int jM, double M, double r);
    int addClumps(cosmology &C, int jz, int jM, double zl, double M, double Sigmac, double r, double phi, rgen &mt, double &kappa, double &gamma1, double &gamma2);
    
};
