#include<iostream>
#include<fstream>
#include<stdlib.h>

#include<adi_solver.h>
#include<observer.h>

using namespace std;


// simulation parameters: to be compared with Parameters1.m in EP11

double LX=0.005;
double LY=0.005;

const int NX=256;
const int NY=256;
const int NS=512;

double deltaZ  =0.001;
double K0 = 2.0*M_PI/1.0e-06;


// for initial condition
comple WhiteNoise(double rx, double ry){ return( comple( rand()/((double) RAND_MAX) -0.5, rand()/((double) RAND_MAX) -0.5 )  );}

int main(int argc, char *argv[]) {

  ofstream fout;
  fout.open("record_space_time.dat");

  ADI_BPM_Solver adi;
  adi.Allocate(LX,NX,LY,NY);
  adi.SetAmplitude(WhiteNoise);
  adi.SetStep(deltaZ, K0);

  // initial condition
  Save1D("sol_i.dat",adi.AM(),adi.CX(),1,adi.NX());

  for(int s=0;s<NS;s++) {


    for(int x=0;x<NX;x++) {
      double ave = 0.0;
      for(int y=0;y<NY;y++) ave += real( adi.AM()[x + y*NX] );
      fout << ave <<" ";
    }
    fout << endl;

    adi.Step();
    cout << s << endl;
  }


  return(0);
}

