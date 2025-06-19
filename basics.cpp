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

// linear interpolation
double lin(double y1, double y2, double x1, double x2, double x) {
    return y1 + (x-x1)/(x2-x1)*(y2 - y1);
}

vector<double> interpolate2(double x0, double y0, vector<double> &x, vector<double> &y, vector<vector<vector<double> > > &f) {
    
    int nx = x.size();
    int ny = y.size();
    int nf = f[0][0].size();
    
    // find the closest corner where x>x0 and y>y0
    int jx = lower_bound(x.begin(), x.end(), x0) - x.begin();
    int jy = lower_bound(y.begin(), y.end(), y0) - y.begin();

    vector<double> tmp(nf, 0.0);
    if (jx <= 0) {
        if (jy <= 0) {
            tmp = f[0][0];
        } else if (jy >= ny) {
            tmp = f[0][ny-1];
        } else {
            for (int j = 0; j < nf; j++) {
                tmp[j] = lin(f[0][jy-1][j], f[0][jy][j], y[jy-1], y[jy], y0);
            }
        }
    } else if (jx >= nx) {
        if (jy <= 0) {
            tmp = f[nx-1][0];
        } else if (jy >= ny) {
            tmp = f[nx-1][ny-1];
        } else {
            for (int j = 0; j < nf; j++) {
                tmp[j] = lin(f[nx-1][jy-1][j], f[nx-1][jy][j], y[jy-1], y[jy], y0);
            }
        }
    } else {
        if (jy <= 0) {
            for (int j = 0; j < nf; j++) {
                tmp[j] = lin(f[jx-1][0][j], f[jx][0][j], x[jx-1], x[jx], x0);
            }
        } else if (jy >= ny) {
            for (int j = 0; j < nf; j++) {
                tmp[j] = lin(f[jx-1][ny-1][j], f[jx][ny-1][j], x[jx-1], x[jx], x0);
            }
        }
        else {
            double low, up;
            for (int j = 0; j < nf; j++) {
                low = lin(f[jx-1][jy-1][j], f[jx][jy-1][j], x[jx-1], x[jx], x0);
                up = lin(f[jx-1][jy][j], f[jx][jy][j], x[jx-1], x[jx], x0);
                tmp[j] = lin(low, up, y[jy-1], y[jy], y0);
            }
        }
    }
    return tmp;
}

void writeToFile(vector<double> &row, const string &filename) {
    ofstream outFile(filename);
    for (int j = 0; j < row.size(); j++) {
        outFile << row[j];
        if (j < row.size() - 1) {
            outFile << endl;
        }
    }
    outFile.close();
}

void writeToFile(vector<vector<double> > &matrix, const string &filename) {
    ofstream outFile(filename);
    for (auto& row : matrix) {
        for (int j = 0; j < row.size(); j++) {
            outFile << row[j];
            if (j < row.size() - 1) {
                outFile << "   ";
            }
        }
        outFile << endl;
    }
    outFile.close();
}

void writeToFile(vector<double> &x, vector<vector<vector<double> > > &cubic, const string &filename) {
    ofstream outFile(filename);
    vector<vector<double> > matrix;
    for (int jx = 0; jx < x.size(); jx++) {
        matrix = cubic[jx];
        for (auto& row : matrix) {
            outFile << x[jx] << "   ";
            for (int j = 0; j < row.size(); j++) {
                outFile << row[j];
                if (j < row.size() - 1) {
                    outFile << "   ";
                }
            }
            outFile << endl;
        }
    }
    outFile.close();
}

void writeToFile(vector<double> &x, vector<double> &y, vector<vector<vector<double> > > &f, const string &filename) {
    ofstream outFile(filename);
    for (int jx = 0; jx < x.size(); jx++) {
        for (int jy = 0; jy < y.size(); jy++) {
            outFile << x[jx] << "   " << y[jy] << "   ";
            for (int j = 0; j < f[jx][jy].size(); j++) {
                outFile << f[jx][jy][j];
                outFile << "   ";
            }
            outFile << endl;
        }
    }
    outFile.close();
}

void writeToFile(vector<double> &x, vector<double> &y, vector<double> &z, vector<vector<vector<vector<double> > > > &f, const string &filename) {
    ofstream outFile(filename);
    for (int jx = 0; jx < x.size(); jx++) {
        for (int jy = 0; jy < y.size(); jy++) {
            for (int jz = 0; jz < z.size(); jz++) {
                outFile << x[jx] << "   " << y[jy] << "   " << z[jz] << "   ";
                for (int j = 0; j < f[jx][jy][jz].size(); j++) {
                    outFile << f[jx][jy][jz][j];
                    outFile << "   ";
                }
                outFile << endl;
            }
        }
    }
    outFile.close();
}

