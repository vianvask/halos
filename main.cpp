#include "cosmology.h"
#include "lensing.h"
#include "UVluminosity.h"

// units: masses in solar masses, time in Myr, length in kpc

int main (int argc, char *argv[]) {
    cout << setprecision(3) << fixed;
    
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
    cout << "t_0 =  " << C.age(0.0) << " Myr." << endl;
    cout << "D_L(z=1) = " << C.DL(1.0) << " kpc." << endl;
    
    double z, M, hmf, dotM, DdotM, MUV, AUV, PhiUV;
    string filename;
    ofstream outfile;
    
    cout << "Computing UV luminosity functions..." << endl;
    
    vector<vector<vector<double> > > UVLFlist = UVLFlistf(C, 1.0e9, 3.2e11, 0.1, -1.0, 0.6);
    
    // output the halo mass function and the UV luminosity function
    filename = "hmf.dat";
    outfile.open(filename.c_str());
    
    for (int jz = 0; jz < C.Nz; jz++) {
        for (int jM = 0; jM < C.NM; jM++) {
            z = C.zlist[jz];
            M = C.hmflist[jz][jM][1];
            hmf = C.hmflist[jz][jM][2];
            dotM = C.dotMlist[jz][jM][2];
            DdotM = C.dotMlist[jz][jM][3];
            
            MUV = UVLFlist[jz][jM][1];
            AUV = UVLFlist[jz][jM][2];
            PhiUV = UVLFlist[jz][jM][3];
            outfile << z << "   " << M << "   " << max(1.0e-64,hmf) << "   " << dotM << "   " << DdotM << "   " << MUV << "   " << AUV << "   " << max(1.0e-64,PhiUV) << endl;
        }
    }
    outfile.close();
    
    vector<double> Zlist {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    vector<vector<double> > Plnmu;
    vector<vector<vector<double> > > Plnmuz;
    
    // read or generate lensing amplification distribution
    ifstream infile;
    infile.open("Plnmu.dat");
    
    if (infile) {
        cout << "Reading lensing amplifications..." << endl;
        
        Zlist.clear();
        vector<double> tmp(2,0.0);
        double A;
        int jA = 0;
        z = 0;
        while (infile >> A) {
            if (jA == 0) {
                if (A > z) {
                    if (z > 0) {
                        Plnmuz.push_back(Plnmu);
                        Plnmu.clear();
                    }
                    z = A;
                    Zlist.push_back(z);
                }
                jA++;
            } else if (jA == 1) {
                tmp[jA-1] = A;
                jA++;
            } else if (jA == 2) {
                tmp[jA-1] = A;
                jA = 0;
            }
            if (jA == 0) {
                Plnmu.push_back(tmp);
            }
        }
        Plnmuz.push_back(Plnmu);
        Plnmu.clear();
    }
    else {
        cout << "Generating lensing amplifications..." << endl;
        
        filename = "Plnmu.dat";
        outfile.open(filename.c_str());
        
        double rS = 10.0; // source radius in kpc
        int Nkappa = 40;
        int Nreal = 20000000;
        int Nbins = 2000;
        
        rgen mt(time(NULL)); // random number generator
        
        for (int jz = 0; jz < Zlist.size(); jz++) {
            z = Zlist[jz];
            cout << "z = " << z << endl;
            Plnmu = Plnmuf(C, Nkappa, z, 1.0, rS, Nreal, Nbins, mt);
            
            // output dP/dlnmu
            for (int jb = 0; jb < Nbins; jb++) {
                outfile << z << "   " << Plnmu[jb][0] << "   " << Plnmu[jb][1] << endl;
            }
            Plnmuz.push_back(Plnmu);
        }
        outfile.close();
    }
    
    cout << "z = { ";
    for (int jz = 0; jz < Zlist.size(); jz++) {
        cout << Zlist[jz] << " ";
    }
    cout << "}" << endl;
    
    cout << "Computing lensed UV luminosity functions..." << endl;
    
    filename = "UVluminosity.dat";
    outfile.open(filename.c_str());
    
    vector<vector<vector<double> > > PhiUVlensed(Zlist.size(), vector<vector<double> > (C.NM, vector<double> (2,0.0)));
    vector<vector<double> > PhiUVlist(C.NM, vector<double> (2,0.0));
    int jz;
    double dlnmu;
    for (int jZ = 0; jZ < Zlist.size(); jZ++) {
        z = Zlist[jZ];
        Plnmu = Plnmuz[jZ];
        
        // find index of z in the longer list of z values
        jz = lower_bound(C.zlist.begin(), C.zlist.end(), z)- C.zlist.begin();
        if (jz > 0 && C.zlist[jz]-z > z-C.zlist[jz-1]) {
            jz--;
        }
                
        // make a list Phi(M_UV) that can be easily interpolated
        for (int jM = 0; jM < C.NM; jM++) {
            PhiUVlist[jM][0] = UVLFlist[jz][jM][1];
            PhiUVlist[jM][1] = UVLFlist[jz][jM][3];
        }
        
        // compute the lensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            MUV = UVLFlist[jz][jM][1];
            AUV = UVLFlist[jz][jM][2];
            
            // convolution integral
            PhiUV = 0.0;
            dlnmu = Plnmu[1][0]-Plnmu[0][0];
            for (int jb = 0; jb < Plnmu.size(); jb++) {
                PhiUV += interpolate(MUV - AUV + 1.08574*Plnmu[jb][0], PhiUVlist)*Plnmu[jb][1]*dlnmu/exp(Plnmu[jb][0]);
            }
            
            PhiUVlensed[jZ][jM][0] = MUV;
            PhiUVlensed[jZ][jM][1] = PhiUV; // dust + lensing
            
            // output the UV luminosity function (no dust + no lensing, dust + no lensing, dust + lensing)
            outfile << z << "   " << MUV << "   " << max(1.0e-64,UVLFlist[jz][jM][3]) << "   " << max(1.0e-64,interpolate(MUV - AUV, PhiUVlist)) << "   " << max(1.0e-64,PhiUV) << endl;
        }
    }
    outfile.close();
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
