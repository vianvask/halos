#include "lensing.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    cout << setprecision(13) << fixed;
    
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
    }
    // run with PDG values
    else {
        C.OmegaM = 0.315;
        C.OmegaB = 0.0493;
        C.zeq = 3402.0;
        C.sigma8 = 0.811;
        C.h = 0.674;
        C.T0 = 2.7255;
        C.ns = 0.965;
    }
    
    // accuracy parameters
    C.Nk = 1000;
    C.NM = 1000;
    C.Nz = 1000;
    
    // mass and redshift ranges
    C.Mmin = 1.0e5;
    C.Mmax = 1.0e17;
    C.zmin = 0.01;
    C.zmax = 20.0;
    
    C.initialize();
    time_req = clock() - time_req;
    cout << "Initialization evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    cout << "The age of the universe today is " << C.age(0.0) << " Myr." << endl;
    cout << "The luminosity distance of a source at z = 1 is " << C.DL(1.0) << " kpc." << endl;
            
    // output the halo mass function
    string filename = "hmf.dat";;
    ofstream outfile;
    outfile.open(filename.c_str());
    double z, M, hmf;
    for (int jz = 0; jz < C.Nz; jz++) {
        for (int jM = 0; jM < C.NM; jM++) {
            z = C.hmflist[jz][jM][0];
            M = C.hmflist[jz][jM][1];
            hmf = C.hmflist[jz][jM][2];
            if (hmf < 1e-64) {
                hmf = 0.0;
            }
            outfile << z << "   " << M << "   " << hmf << endl;
        }
    }
    outfile.close();
    
    vector<vector<double> > N1list = dNdlnmu(C, 20, 2.0, 0.001, 1.0);
    
    // output the distribution dN_1/dlnmu
    filename = "N1.dat";;
    outfile.open(filename.c_str());
    outfile << setprecision(12) << fixed;
    for (int j = 0; j < N1list.size(); j++) {
        outfile << N1list[j][0] << "   " << N1list[j][1] << "   " << N1list[j][2] << endl;
    }
    outfile.close();
        
    double Nh = N1list[N1list.size()-1][2];
    cout << "Number of halos along each path: " << Nh << endl;

    vector<vector<double> > Fm1(N1list.size(), vector<double> (2,0.0));
    for (int j = 0; j < Fm1.size(); j++) {
        Fm1[j][0] = N1list[j][2]/Nh;
        Fm1[j][1] = log(N1list[j][0]);
    }
    
    // random number generator
    rgen mt(time(NULL));
    
    double x, lnmutot;
    vector<double> lnmu(1000000, 0.0);
    filename = "lnmu.dat";
    outfile.open(filename.c_str());
    for (int j = 0; j<lnmu.size(); j++) {
        lnmutot = 0.0;
        for (int jh = 0; jh < Nh; jh++) {
            x = randomreal(0.0,1.0,mt);
            lnmutot += exp(interpolate(x,Fm1));
        }
        outfile << lnmutot << endl;
    }
    outfile.close();
    
    time_req = clock() - time_req;
    cout << "Evaluation time after initialization: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
