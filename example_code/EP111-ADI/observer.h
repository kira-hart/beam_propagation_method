#ifndef __observer_h__
#define __observer_h__

#include<iostream>
#include<fstream>
#include<string>
#include<complex>

using namespace std;

typedef complex<double> comple;

void Save2D(string filename, const comple *src,  int nx, int ny, double (*function)(comple));

void Save1D(string filename, const comple *src, const double *c, int stride, int n);

double fun_real(comple x);
double fun_imag(comple x);
double fun_norm(comple x);
double fun_phas(comple x);

#endif // __observer_h__
