#include<iostream>
#include<fstream>
#include<stdlib.h>

#include<adi_solver.h>
#include<observer.h>

using namespace std;

// choose parameter file and re-compile
#include"parameters1.cc"

// for initial condition
comple gaussxy(double wr, double x, double y, double z, double k, double f, double kx, double ky);
comple GaussBeam(double rx, double ry){ return( gaussxy( W0,  rx,  ry,  currentZ,  K0,  F, 0.0,  0.0) );}

int main(int argc, char *argv[]) {

  ADI_BPM_Solver adi;
  adi.Allocate(LX,NX,LY,NY);
  adi.SetAmplitude(GaussBeam);
  adi.SetStep(deltaZ, K0);

  // initial condition
  Save1D("sol_i.dat",&adi.AM()[NX*NY/2],adi.CX(),1,adi.NX());

  for(;currentZ<finalZ;) {
    adi.Step();
    currentZ += deltaZ;
    //    cout <<"curentZ: " << currentZ << endl;
  }

  // record final solution data
  Save1D("sol_f.dat",&(adi.AM()[NY*NX/2]),adi.CX(),1,adi.NX());
  Save2D("map_f.dat",  adi.AM(), NX, NY, fun_real);

  // for error evaluation:
  comple err[NX*NY];
  memcpy(err,adi.AM(),NX*NY*sizeof(comple));

  // analytic solution for comparison
  adi.SetAmplitude(GaussBeam);
  Save1D("sol_a.dat",&(adi.AM()[NY*NX/2]),adi.CX(),1,adi.NX());
  Save2D("map_a.dat",  adi.AM(), NX, NY, fun_real);

  // subtract numerical and target 
  for(int i=0;i<NX*NY;i++) err[i] -= adi.AM()[i];
  Save1D("sol_e.dat",&(err[NY*NX/2]),adi.CX(),1,adi.NX());
  Save2D("map_e.dat",  err, NX, NY, fun_norm);

  return(0);
}


comple gaussxy(double wr, double x, double y, double z, double k, double f, double kx, double ky) 
{
  comple aux1, aux2, imi(0.0,1.0);

  if( wr <= 0.0 ) {
    // constant profile
    return(1.0);
  }

  if(f != 0.0) 
    aux1 =  1.0/(wr*wr) + imi*k/(2.0*f);
  else
    aux1 =  1.0/(wr*wr);

  aux2 =  1.0/(1.0 + 2.0*imi*z/k*aux1);

  return( aux2*exp(-(x*x+y*y)*aux1*aux2)*exp(imi*(x*kx + y*ky)));
}


