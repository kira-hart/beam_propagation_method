#ifndef __bpm_icn_h__
#define __bpm_icn_h__

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<complex>
#include"uppe_constants.h"
#include"bpm_parameters.h"


using namespace std;

typedef complex<double> comple;

#define index(x,y) ((x) + NX*(y))


class ICN {

public:

  Params_All params;

  int NX;    // transverse dimension
  int NY;    // transverse dimension
  int NC;    // num of  polarization components (usually=2) 
  int NA;    // num of auxiliay fields (to implement mid-steps)

  comple *data;       // pointer to array holding all beam amplitudes
  comple *pmls;       // PMLayers
  comple *pml_aa_x;   // PML coefficient (contained in pmls)
  comple *pml_bb_x;
  comple *pml_aa_y;
  comple *pml_bb_y;

  double *epsr;       // realtive permitivity map
  double *cx;         // coordinates
  double *cy;         // coordinates

  double dz;          // integration step
  double dx;          // grid spacings
  double dy;

  double k0;          // free-space wavenumber
  double nref;        // referencde refractive index

  // coefficients in discretized equations:
  comple coeff_xx;    // laplacian along x,...
  comple coeff_yy;    // along y
  comple coeff_xy;    // mixed derivatives 
  comple coeff_di;    // this is the diagonal

  ICN();
 ~ICN();

 // initialization
 void Allocate(int nx, int ny, double LX, double LY);
 void SetStep(comple delta_z, double lambda);
 void Normalize();

 // pml
 void PMLinit(bool on);

 // access to fields, physical and auxiliary
 comple* EX(int c) { return( &data[NX*NY*NC*c        ] ); }
 comple* EY(int c) { return( &data[NX*NY*NC*c + NX*NY] ); }

 // function implemeting discretized propagation equation
 void ApplyAXX(int src, int tar);
 void ApplyAYY(int src, int tar);
 void ApplyAXY(int src, int tar);
 void ApplyAYX(int src, int tar);
 void ApplyRHS(int src, int tar);

 // stepping 
 void Add(int src, int tar, double coeff);
 void Step();
 
 // observer functions
 void Report(int c, int n);
 void SaveMatrix(const char *fname, int part, comple *src);

 // geomtery
 double EpsilonRelative(double rx, double ry);
};

#endif // __bpm_icn_h__
