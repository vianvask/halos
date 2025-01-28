#include "lensing.h"

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
    
    vector<vector<vector<double> > > UVLFlist = C.UVLFlistf(1.0e9, 3.2e11, 0.1, -1.2, 0.5);
            
    // output the halo mass function and the UV luminosity function
    string filename = "hmf.dat";;
    ofstream outfile;
    outfile.open(filename.c_str());
    double z, M, hmf, dotM, DdotM, MUV, AUV, phiUV;
    for (int jz = 0; jz < C.Nz; jz++) {
        for (int jM = 0; jM < C.NM; jM++) {
            z = C.hmflist[jz][jM][0];
            M = C.hmflist[jz][jM][1];
            hmf = C.hmflist[jz][jM][2];
            dotM = C.dotMlist[jz][jM][2];
            DdotM = C.dotMlist[jz][jM][3];
            MUV = UVLFlist[jz][jM][2];
            AUV = UVLFlist[jz][jM][3];
            phiUV = UVLFlist[jz][jM][4];
            if (hmf < 1e-64) {
                hmf = 0.0;
            }
            if (abs(phiUV) < 1e-64) {
                phiUV = 0.0;
            }
            outfile << z << "   " << M << "   " << hmf << "   " << dotM << "   " << DdotM << "   " << MUV << "   " << AUV << "   " << phiUV << endl;
        }
    }
    outfile.close();
    
    cout << "Generating lensing amplifications..." << endl;
            
    int Nkappa = 80, Nreal = 100000000, Nbins = 10000;
    double Nh, x, kappatot, lnmutot, lnmumin, lnmumax, lnmumean, dlnmu, dP = 1.0/(1.0*Nreal);
    
    vector<double> zlist {0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0};

    vector<vector<double> > N1list;
    vector<vector<double> > Fm1(Nkappa, vector<double> (2,0.0));
    vector<double> lnmu(Nreal, 0.0);
    vector<vector<double> > Plnmu(Nbins, vector<double> (2, 0.0));
    
    filename = "Plnmu.dat";
    outfile.open(filename.c_str());
    
    // random number generator
    rgen mt(time(NULL));
    
    for (int jz = 0; jz < zlist.size(); jz++) {
        cout << "z = " << zlist[jz] << endl;
        
        // source size = 10.0 kpc
        N1list = dNdkappa(C, Nkappa, zlist[jz], 1.0, 10.0);
        Nh = N1list[N1list.size()-1][2];
        
        for (int j = 0; j < Fm1.size(); j++) {
            Fm1[j][0] = N1list[j][2]/Nh;
            Fm1[j][1] = log(N1list[j][0]);
        }
        
        // generate realizations of kappa
        lnmumin = 1.0;
        lnmumax = 0.0;
        lnmumean = 0.0;
        for (int j = 0; j<Nreal; j++) {
            kappatot = 0.0;
            for (int jh = 0; jh < Nh; jh++) {
                x = randomreal(0.0,1.0,mt);
                kappatot += exp(interpolate(x,Fm1));
            }
            
            lnmutot = log(pow(1.0-kappatot,-2.0)); // mu = (1-kappa)^-2 (see 1106.3823)
            if (lnmutot > lnmumax) {
                lnmumax = lnmutot;
            }
            if (lnmutot < lnmumin) {
                lnmumin = lnmutot;
            }
            lnmumean += lnmutot;
            lnmu[j] = lnmutot;
        }
        lnmumean = lnmumean/(1.0*Nreal);
        
        // binning
        dlnmu = (lnmumax-lnmumin)/(1.0*(Nbins-1));
        for (int j = 0; j < Nbins; j++) {
            Plnmu[j][0] = (lnmumin-lnmumean) + j*dlnmu;
        }
        for (int j = 0; j < Nreal; j++) {
            Plnmu[(int) round((Nbins-1)*(lnmu[j]-lnmumin)/(lnmumax-lnmumin))][1] += dP/dlnmu;
        }
        
        // output P(kappa)
        for (int j = 0; j < Nbins; j++) {
            outfile << zlist[jz] << "   " << Plnmu[j][0] << "   " << Plnmu[j][1] << endl;
            Plnmu[j][1] = 0.0;
        }
    }
    outfile.close();
    
    time_req = clock() - time_req;
    cout << "Total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " min." << endl;
    
    return 0;
}
