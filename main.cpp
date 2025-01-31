#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    cout << setprecision(3) << fixed;
    
    cout << "Preparing..." << endl;
    
    // timing
    clock_t time_req = clock();

    // input parameters: PDG values
    cosmology C;
    C.OmegaM = 0.315;
    C.OmegaB = 0.0493;
    C.zeq = 3402.0;
    C.sigma8 = 0.811;
    C.h = 0.674;
    C.T0 = 2.7255;
    C.ns = 0.965;
    
    // accuracy parameters
    C.Nk = 1000;
    C.NM = 1200;
    C.Nz = 100;
    
    // mass and redshift ranges
    C.Mmin = 1.0e5;
    C.Mmax = 1.0e17;
    C.zmin = 0.01;
    C.zmax = 20.01;
    
    C.initialize();
    //cout << "t_0 =  " << C.age(0.0) << " Myr." << endl;
    //cout << "D_L(z=1) = " << C.DL(1.0) << " kpc." << endl;
    
    double z, M, hmf, dotM, DdotM;
    ofstream outfile;
    
    // output the halo mass function and the UV luminosity function
    outfile.open("hmf.dat");
    for (int jz = 0; jz < C.Nz; jz++) {
        for (int jM = 0; jM < C.NM; jM++) {
            z = C.zlist[jz];
            M = C.hmflist[jz][jM][1];
            hmf = C.hmflist[jz][jM][2];
            dotM = C.dotMlist[jz][jM][2];
            DdotM = C.dotMlist[jz][jM][3];
            
            outfile << z << "   " << M << "   " << max(1.0e-64,hmf) << "   " << dotM << "   " << DdotM << endl;
        }
    }
    outfile.close();
    
    // list of z values for lensing
    vector<double> Zlist {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    
    // lensing source radius in kpc
    double rS = 10.0;
    
    // accuracy parameters for lensing
    int Nkappa = 40;
    int Nreal = 20000000;
    int Nbins = 2000;
    
    // read or generate lensing amplification distribution
    vector<vector<vector<double> > > Plnmuz = getPlnmu(C, Zlist, rS, Nkappa, Nreal, Nbins);
    
    cout << "Computing UV luminosity functions..." << endl;
    
    // benchmark case
    vector<vector<vector<double> > > PhiUVlensed = PhiUV(C, Zlist, Plnmuz, 1.0e9, 1.4e11, 0.06, -1.0, 0.15, 1.0);
    
    // output the UV luminosity function (no dust + no lensing, dust + no lensing, dust + lensing)
    outfile.open("UVluminosity.dat");
    for (int jZ = 0; jZ < Zlist.size(); jZ++) {
        z = Zlist[jZ];
        for (int jM = 0; jM < C.NM; jM++) {
            outfile << z << "   " << PhiUVlensed[jZ][jM][0] << "   " << max(1.0e-64,PhiUVlensed[jZ][jM][1]) << "   " << max(1.0e-64,PhiUVlensed[jZ][jM][2]) << "   " << max(1.0e-64,PhiUVlensed[jZ][jM][3]) << endl;
        }
    }
    outfile.close();
    
    // read UV luminosity data files
    vector<string> datafiles {"UVLF_ 2102.07775.txt", "UVLF_ 2108.01090.txt"};
    vector<vector<double> > data = readUVdata(datafiles);
    cout << "Read " << data.size() << " UV luminosity data points." << endl;
    
    cout << "Computing UV luminosity fit..." << endl;
    
    // TODO: scan over the SFR parameters
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
