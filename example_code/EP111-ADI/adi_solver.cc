#include<iostream>
#include<adi_solver.h>
#include<tdm_utils.h>

using namespace std;

void ADI_BPM_Solver::Allocate(double LX, int NX, double LY, int NY) {
   // this and other functions simplified: no checking and other important stuff

   // define corrdinates
   dimX = NX;
   dimY = NY;
   cx   = new double[dimX];
   cy   = new double[dimY];
   for(int i=0;i<dimX;i++) cx[i] = (LX/dimX)*(i - dimX/2.0);
   for(int i=0;i<dimY;i++) cy[i] = (LY/dimY)*(i - dimY/2.0);

   // 2D array holders
   const int num2D=6;
   Holder2D = new comple[dimX*dimY*num2D];
   E00 = &Holder2D[0*dimX*dimY];
   E12 = &Holder2D[1*dimX*dimY];
   E11 = &Holder2D[2*dimX*dimY];
   R00 = &Holder2D[3*dimX*dimY];
   R12 = &Holder2D[4*dimX*dimY];

   // 1D array holders
   const int num1D=7;
   Holder1DX = new comple[dimX*num1D];
   Holder1DY = new comple[dimY*num1D];

   LPAX = &Holder1DX[0*dimX];
   LPBX = &Holder1DX[1*dimX];
   LPCX = &Holder1DX[2*dimX];
   LMAX = &Holder1DX[3*dimX];
   LMBX = &Holder1DX[4*dimX];
   LMCX = &Holder1DX[5*dimX];
   auxX = &Holder1DX[6*dimX];


   LPAY = &Holder1DY[0*dimY];
   LPBY = &Holder1DY[1*dimY];
   LPCY = &Holder1DY[2*dimY];
   LMAY = &Holder1DY[3*dimY];
   LMBY = &Holder1DY[4*dimY];
   LMCY = &Holder1DY[5*dimY];
   auxY = &Holder1DY[6*dimY];

 }

void ADI_BPM_Solver::SetStep(double dz, double k0) {
  // assuming homogeneous (equidistant) grid(s)
  double dx = cx[1]-cx[0];
  double dy = cy[1]-cy[0];

  comple ideltaX;
  comple ideltaY;
  ideltaX = comple(0.0,1.0)*dz/(4.0*k0*dx*dx);
  ideltaY = comple(0.0,1.0)*dz/(4.0*k0*dy*dy);

  int i;
  for(i=0;i<dimX;i++) LPAX[i] =     +1.0*ideltaX; 
  for(i=0;i<dimX;i++) LPBX[i] = 1.0 -2.0*ideltaX; 
  for(i=0;i<dimX;i++) LPCX[i] =     +1.0*ideltaX; 

  for(i=0;i<dimX;i++) LMAX[i] =     -1.0*ideltaX; 
  for(i=0;i<dimX;i++) LMBX[i] = 1.0 +2.0*ideltaX; 
  for(i=0;i<dimX;i++) LMCX[i] =     -1.0*ideltaX; 

  for(i=0;i<dimX;i++) auxX[i] =  LMBX[i]; 

  for(i=0;i<dimY;i++) LPAY[i] =     +1.0*ideltaY; 
  for(i=0;i<dimY;i++) LPBY[i] = 1.0 -2.0*ideltaY; 
  for(i=0;i<dimY;i++) LPCY[i] =     +1.0*ideltaY; 

  for(i=0;i<dimY;i++) LMAY[i] =     -1.0*ideltaY; 
  for(i=0;i<dimY;i++) LMBY[i] = 1.0 +2.0*ideltaY; 
  for(i=0;i<dimY;i++) LMCY[i] =     -1.0*ideltaY; 

  for(i=0;i<dimY;i++) auxY[i] =  LMBY[i]; 
}


void ADI_BPM_Solver::SetAmplitude(comple (*func) (double rx, double ry)) {
  int x,y;
  for(y=0;y<dimY;y++) for(x=0;x<dimX;x++) {
      E00[x+y*dimX] = func(cx[x],cy[y]);
    }
}


void ADI_BPM_Solver::Step() {
  int x,y;

  // R00  =  LPX*E00
  for(y=0;y<dimY;y++) {
    TDM_Apply(LPAX,LPBX,LPCX, &E00[0 + y*dimX], &R00[0 + y*dimX],    1, dimX);
  }
  
  // Solve:  LMY*E12 = R00,  R00 destroyed
  for(x=0;x<dimX;x++) {
    TDM_Solve(LMAY,LMBY,LMCY, &R00[x + 0*dimX], &E12[x + 0*dimX], dimX, dimY);
    // restore diagonal that was destroyed
    for(y=0;y<dimY;y++) LMBY[y] = auxY[y];
    memcpy(LMBY,auxY,dimY*sizeof(comple));
  }
  
  // R12[x + 0*dimX]
  for(x=0;x<dimX;x++) {
    TDM_Apply(LPAY,LPBY,LPCY, &E12[x + 0*dimX], &R12[x + 0*dimX], dimX, dimY);
  }

  // Solve:  LMX*E11 = R12, R12 destroyed
  for(y=0;y<dimY;y++) {
    TDM_Solve(LMAX,LMBX,LMCX, &R12[0 + y*dimX], &E11[0 + y*dimX],    1, dimX);
    // restore diagonal that was destroyed
    memcpy(LMBX,auxX,dimX*sizeof(comple));
  }
  
  memcpy(E00,E11,dimX*dimY*sizeof(comple));
}
