#include <complex>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

typedef mt19937_64 rgen;

const double PI = 3.141592653589793238463;

string to_string_prec(const double a, const int n);

double randomreal(double x1, double x2, rgen &mt);
double randomreal(double x1, double x2);

double interpolate(double X, vector<vector<double> > &xy);
vector<double> interpolaten(double X, vector<vector<double> > &xy);
vector<double> interpolate2(double x0, double y0, vector<double> &x, vector<double> &y, vector<vector<vector<double> > > &f);

void writeToFile(vector<double> &row, const string &filename);
void writeToFile(vector<vector<double> > &matrix, const string &filename);
void writeToFile(vector<double> &x, vector<vector<vector<double> > > &cubic, const string &filename);
void writeToFile(vector<double> &x, vector<double> &y, vector<vector<vector<double> > > &f, const string &filename);
void writeToFile(vector<double> &x, vector<double> &y, vector<double> &z, vector<vector<vector<vector<double> > > > &f, const string &filename);
