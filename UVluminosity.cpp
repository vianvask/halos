#include "cosmology.h"
#include "UVluminosity.h"

double kappaUV = 1.15e-22; // Msun s erg^-1 Myr^-1

// star formation rate
double fstar(double M, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return epsilon*pow(pow(M/Mc,alpha/gamma) + pow(M/Mc,beta/gamma),-gamma);
}

// derivative of the star formation rate, df_*/dM
double Dfstar(double M, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return -(alpha/(1+pow(M/Mc,(-alpha+beta)/gamma)) + beta/(1+pow(M/Mc,(alpha-beta)/gamma)))*fstar(M,Mc,epsilon,alpha,beta,gamma)/M;
}
                       
                       

// the UV magnitude
double MUV(double M, double dotM, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return 51.63 - 1.08574*log(fstar(M,Mc,epsilon,alpha,beta,gamma)*dotM/kappaUV);
}
    
// derivative of the UV magnitude, dM_UV/dM
double DMUV(double M, double dotM, double DdotM, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return -1.08574*(DdotM/dotM + Dfstar(M,Mc,epsilon,alpha,beta,gamma)/fstar(M,Mc,epsilon,alpha,beta,gamma));
}

// dust extinction, MUV -> MUV - AUV (see 1406.1503 and Table 3 of 1306.2950)
double AUV(double MUV, double z) {
    double C0 = 4.4, C1 = 2.0, sigmabeta = 0.34;
    double beta = -exp(0.17*(19.5+MUV)/(1.54+0.075*z))*(1.54+0.075*z);
    return max(0.0, C0 + 0.2*log(10)*pow(C1*sigmabeta,2.0) + C1*beta);
}

// UV luminosity function as a function of halo mass M (see e.g. 1906.06296)
double UVLF(double M, double dotM, double DdotM, double dndlnM, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma) {
    return -exp(-Mt/M)*dndlnM/M/DMUV(M,dotM,DdotM,Mc,epsilon,alpha,beta,gamma);
}


// UV luminosity function, {z, M_UV, Phi(no dust + no lensing), Phi(dust + no lensing), Phi(dust + lensing)}
vector<vector<vector<double> > > PhiUV(cosmology &C, vector<double> &Zlist, vector<vector<vector<double> > > &Plnmuz, double Mt, double Mc, double epsilon, double alpha, double beta, double gamma) {
    
    vector<vector<vector<double> > > PhiUVlensed(Zlist.size(), vector<vector<double> > (C.NM, vector<double> (4,0.0)));
    vector<vector<double> > PhiUVlist(C.NM, vector<double> (3,0.0));
    
    int jzp;
    double z, dlnmu, Mj, dotMj, DdotMj, dndlnMj, MUVj, AUVj, PhiUVj;
    vector<vector<double> > Plnmu(Plnmuz[0].size(), vector<double> (2,0.0));
    for (int jZ = 0; jZ < Zlist.size(); jZ++) {
        z = Zlist[jZ];
        Plnmu = Plnmuz[jZ];
        
        // find index of z in the longer list of z values
        jzp = lower_bound(C.zlist.begin(), C.zlist.end(), z)- C.zlist.begin();
        if (jzp > 0 && C.zlist[jzp]-z > z-C.zlist[jzp-1]) {
            jzp--;
        }
             
        // compute the unlensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            Mj = C.dotMlist[jzp][jM][1];
            dotMj = C.dotMlist[jzp][jM][2];
            DdotMj = C.dotMlist[jzp][jM][3];
            dndlnMj = C.hmflist[jzp][jM][2];
            
            MUVj = MUV(Mj, dotMj, Mc, epsilon, alpha, beta, gamma);
            
            PhiUVlist[jM][0] = MUVj;
            PhiUVlist[jM][1] = AUV(MUVj, z);
            PhiUVlist[jM][2] = UVLF(Mj, dotMj, DdotMj, dndlnMj, Mt, Mc, epsilon, alpha, beta, gamma);
        }
        
        // compute the lensed UV luminosity function
        for (int jM = 0; jM < C.NM; jM++) {
            MUVj = PhiUVlist[jM][0];
            AUVj = PhiUVlist[jM][1];
            
            // convolution integral
            PhiUVj = 0.0;
            dlnmu = Plnmu[1][0]-Plnmu[0][0];
            for (int jb = 0; jb < Plnmu.size(); jb++) {
                PhiUVj += interpolaten(MUVj - AUVj + 1.08574*Plnmu[jb][0], PhiUVlist)[2]*Plnmu[jb][1]*dlnmu/exp(Plnmu[jb][0]);
            }
            
            PhiUVlensed[jZ][jM][0] = MUVj;
            PhiUVlensed[jZ][jM][1] = PhiUVlist[jM][2]; // no dust + no lensing
            PhiUVlensed[jZ][jM][2] = interpolaten(MUVj - AUVj, PhiUVlist)[2]; // dust + lensing
            PhiUVlensed[jZ][jM][3] = PhiUVj; // dust + lensing
        }
    }
    return PhiUVlensed;
}

// read UV luminosity data, {MUV, Phi, +sigmaPhi, -sigmaPhi}
vector<vector<double> > readUVdata(string filename) {
    ifstream infile;
    infile.open(filename);
    vector<vector<double> > data;
    vector<double> line(4,0.0);
    int jline = 0;
    double A;
    if (infile) {
        while (infile >> A) {
            line[jline] = A;
            jline++;
            if (jline == 4) {
                jline = 0;
                data.push_back(line);
            }
        }
    } else {
        cout << "couldn't open " << filename << endl;
    }
    infile.close();
    
    return data;
}
vector<vector<double> > readUVdata(vector<string> filenames) {
    vector<vector<double> > data;
    vector<vector<double> > data0;
    
    for (int j = 0; j < filenames.size(); j++) {
        data0 = readUVdata(filenames[j]);
        data.insert(data.end(), data0.begin(), data0.end());
    }
    return data;
}
