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
double randomreal(double x1, double x2) {
    long double r01 = rand()/(1.0*RAND_MAX);
    return (x1 + (x2-x1)*r01);
}


// linear interpolation
// xy[all][0] must be strictly increasing or strictly decreasing
double interpolate(double X, vector<vector<double> > &xy) {
    
    int n = xy.size();
    auto comp0 = [&](vector<double> &XY, double val) {
      return XY[0] < val;
    };
    
    int jx;
    if (xy[0][0] < xy[n-1][0]) {
        jx = lower_bound(xy.begin(), xy.end(), X, comp0) - xy.begin();
    } else {
        jx = xy.rend() - lower_bound(xy.rbegin(), xy.rend(), X, comp0);
    }
    
    double Y;
    if (jx <= 0) {
        Y = xy[0][1];
    } else if (jx >= n) {
        Y = xy[n-1][1];
    } else {
        Y = xy[jx-1][1] + (X-xy[jx-1][0])/(xy[jx][0]-xy[jx-1][0])*(xy[jx][1]-xy[jx-1][1]);
    }
    return Y;
}

// linear interpolation of vector of form {{x_1,y_11,y_21,...},{x_2,y_12,y_22,...},...}
// xy[all][0] must be strictly increasing or strictly decreasing
vector<double> interpolaten(double X, vector<vector<double> > &xy) {
    
    int n = xy.size();
    auto comp0 = [&](vector<double> &XY, double val) {
      return XY[0] < val;
    };
    
    int jx;
    if (xy[0][0] < xy[n-1][0]) {
        jx = lower_bound(xy.begin(), xy.end(), X, comp0) - xy.begin();
    } else {
        jx = xy.rend() - lower_bound(xy.rbegin(), xy.rend(), X, comp0);
    }
    
    vector<double> tmp(xy[0].size(),0.0);
    if (jx <= 0) {
        tmp = xy[0];
    } else if (jx >= n) {
        tmp = xy[n-1];
    } else {
        tmp[0] = X;
        for (int j = 1; j < tmp.size(); j++) {
            tmp[j] = xy[jx-1][j] + (X-xy[jx-1][0])/(xy[jx][0]-xy[jx-1][0])*(xy[jx][j]-xy[jx-1][j]);
        }
    }
    return tmp;
}

