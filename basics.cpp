#include "basics.h"

string to_string_prec(double a, const int n) {
    ostringstream out;
    out.precision(n);
    out << fixed << a;
    return out.str();
}

// random real number in the range (x_1,x_2)
double randomreal(double x1, double x2, rgen &mt) {
    long double r01 = mt()/(1.0*mt.max());
    return (x1 + (x2-x1)*r01);
}

// linear interpolation
double interpolate(double x, vector<vector<double> > &y) {
    int n = y.size();
    if (x > y[n-1][0]) {
        //cout << "Warning: the point lies above of the interpolation range." << endl;
        return y[n-1][1];
    }
    if (x < y[0][0]) {
        //cout << "Warning: the point lies below of the interpolation range." << endl;
        return y[0][1];
    }
    
    int jmin = 0, jmax = n-1, jx;
    while (jmax > jmin+1) {
        jx = (int) floor((jmax+jmin)/2.0);
        if (y[jx][0] < x) {
            jmin = jx;
        } else {
            jmax = jx;
        }
    }
    jx = jmin;
    return y[jx][1] + (x - y[jx][0])/(y[jx+1][0] - y[jx][0])*(y[jx+1][1] - y[jx][1]);
}
vector<double> interpolaten(double x, vector<vector<double> > &y) {
    int n = y.size();
    if (x > y[n-1][0]) {
        //cout << "Warning: the point lies above of the interpolation range." << endl;
        return y[n-1];
    }
    if (x < y[0][0]) {
        //cout << "Warning: the point lies below of the interpolation range." << endl;
        return y[0];
    }
    
    int jmin = 0, jmax = n-1, jx;
    while (jmax > jmin+1) {
        jx = (int) floor((jmax+jmin)/2.0);
        if (y[jx][0] < x) {
            jmin = jx;
        } else {
            jmax = jx;
        }
    }
    jx = jmin;
    
    vector<double> tmp(y[0].size(),0.0);
    tmp[0] = x;
    for (int j = 1; j < tmp.size(); j++) {
        tmp[j] = y[jx][j] + (x - y[jx][0])/(y[jx+1][0] - y[jx][0])*(y[jx+1][j] - y[jx][j]);
    }
    return tmp;
}
