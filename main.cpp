#include "functions.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    cout << setprecision(15) << fixed;
    
    // timing
    clock_t time_req = clock();
    
    // input parameters
    cosmology C;
    if (argc > 1) {
        C.OmegaM = atof(argv[1]);
        C.OmegaB = atof(argv[2]);
        C.zeq = atof(argv[3]);
        C.sigma8 = atof(argv[4]);
        C.h = atof(argv[5]);
        C.T0 = atof(argv[6]);
        C.ns = atof(argv[7]);
    } else { // run with PDG values:
        C.OmegaM = 0.315;
        C.OmegaB = 0.0493;
        C.zeq = 3402.0;
        C.sigma8 = 0.811;
        C.h = 0.674;
        C.T0 = 2.7255;
        C.ns = 0.965;
    }
    
    // derived parameters
    C.OmegaR = C.OmegaM/(1+C.zeq);
    C.OmegaL = 1.0 - C.OmegaM - C.OmegaR;
    C.OmegaC = C.OmegaM - C.OmegaB;
    C.H0 = 0.000102247*C.h;
    C.rhoc = 277.394*pow(C.h,2.0);
    C.rhoM0 = C.OmegaM*C.rhoc;
    
    // accuracy parameters
    int Nk = 1000;
    int NM = 1000;
    int Nz = 1000;
    double Mmin = 1.0, Mmax = 1.0e16;
    double zmin = 0.01, zmax = 20.0;
    
    // fix deltaH to match the input sigma8
    double deltaH = C.sigma8/C.sigma(8000.0/C.h, 1.0, Nk);
    
    // compute the variance of matter fluctuations, {M,sigma(M),sigma'(M)}
    double dlogM = (log(Mmax)-log(Mmin))/(1.0*NM);
    vector<vector<double> > sigma(NM, vector<double> (3,0.0));
    double RM, M = Mmin;
    double sigman = 0.0, sigmanp = sigman;
    for (int jM = 0; jM < NM; jM++) {
        RM = pow(3.0*M/(4.0*PI*C.rhoM0),1.0/3.0);
        
        sigmanp = sigman;
        sigman = C.sigma(RM, deltaH, Nk);
        
        sigma[jM][0] = M;
        sigma[jM][1] = sigman;
        sigma[jM][2] = (sigman-sigmanp)/(exp(log(M)+dlogM)-M);

        M = exp(log(M)+dlogM);
    }
    sigma[0][2] = sigma[1][2];
    
    // compute the Seth-Tormen halo mass function
    function<double(double)> nuf = [&](double nu) {
        double p = 0.3;
        double q = 0.75;
        double A = 1.0/(1+(pow(2.0,-p)*tgammaf(0.5-p)/sqrt(PI)));
        return A*(1+pow(q*nu,-p))*sqrt(q*nu/(2.0*PI))*exp(-q*nu/2.0);
    };
    double dlogz = (log(zmax)-log(zmin))/(1.0*Nz);
    vector<vector<vector<double> > > dndlnM(NM, vector<vector<double> > (Nz, vector<double> (3,0.0)));
    double z = zmin;
    for (int jz = 0; jz < Nz; jz++) {
        for (int jM = 0; jM < NM; jM++) {
            M = sigma[jM][0];
            
            dndlnM[jz][jM][0] = z;
            dndlnM[jz][jM][1] = M;
            dndlnM[jz][jM][2] = -C.rhoM0*nuf(pow(C.deltac(z)/sigma[jM][1],2.0))*2.0*sigma[jM][2]/sigma[jM][1];
        }
        z = exp(log(z) + dlogz);
    }
    
    // output the halo mass function
    string filename = "hmf.dat";;
    ofstream outfile;
    outfile.open(filename.c_str());
    double hmf;
    for (int jz = 0; jz < Nz; jz++) {
        for (int jM = 0; jM < NM; jM++) {
            z = dndlnM[jz][jM][0];
            M = dndlnM[jz][jM][1];
            hmf = dndlnM[jz][jM][2];
            if (hmf < 1e-64) {
                hmf = 0.0;
            }
            outfile << z << "   " << M << "   " << hmf << endl;
        }
    }
    outfile.close();
            
    // random number generator
    rgen mt(time(NULL)*(C.OmegaM));
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
