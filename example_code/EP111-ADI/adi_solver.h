#ifndef __adi_solver_h__
#define __adi_solver_h__

#include<complex>
typedef std::complex<double> comple;

class ADI_BPM_Solver {

 private:

  // grid dimensions
  int dimX;
  int dimY;

  // grid coordinates (eventually will become complex-valued for PML)
  double *cx;
  double *cy;

  // 2D arrays
  comple *Holder2D;
  comple *E00;      // current amplitude
  comple *E12;      // half-step amplitude
  comple *E11;      // full-step, new amplitude
  comple *R00;      // first half-step right-hand-side
  comple *R12;      // second half-step right-hand-side

/*
the ADI scheme:
src -> E00
TDM_Apply:  R00  =  LPX*E00
TDM_Solve:          LMY*E12 = R00
TDM_APPLY:  R12  =  LPY*E12
TDM_Solve:          LMX*E11 = R12
E11-> tar
*/

  // 1D arrays: hold diagonals of C-N matrices
  comple *Holder1DX;
  comple *Holder1DY;

  // matrix L+ along X
  comple *LPAX;   // diagonal a
  comple *LPBX;   // diagonal b
  comple *LPCX;   // diagonal c

  // matrix L- along X
  comple *LMAX;
  comple *LMBX;
  comple *LMCX;

  // matrix L+ along Y
  comple *LPAY;
  comple *LPBY;
  comple *LPCY;

  // matrix L- along Y
  comple *LMAY;
  comple *LMBY;
  comple *LMCY;

  // 1D auxiliaries: source vectors and targets to store destroyed diagonals b
  comple *auxX;
  comple *auxY;

 public:
  
  ADI_BPM_Solver() {cx=cy=0; Holder1DX=Holder1DY=0; Holder2D=0;};
 ~ADI_BPM_Solver() {delete[] cx; delete[] cy; delete[] Holder1DX; delete Holder1DY; delete[] Holder2D;};

  void Allocate(double LX, int NX, double LY, int NY);
  void SetAmplitude(comple (*Function)(double rx, double ry));
  void SetStep(double deltaZ, double k0);
  void Step();  

  // pass handles to observer
  const double* CX() {return(cx); }
  const double* CY() {return(cy); }
  const comple* AM() {return(E00);}
  int NX() {return(dimX);}
  int NY() {return(dimY);}
};

#endif //__adi_solver_h__
