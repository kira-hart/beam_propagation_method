#include<stdio.h>
#include<cstdlib>
#include<complex>

using namespace std;

#include"uppe_solver.h"

const int NX = 300;

complex<double> AA[NX];
complex<double> BB[NX];

class RHSpars {

public:
  int    nx;
  double la;
  double nr;
  int    iiL;
  int    iiR;
  double dr;
};

void SetInterfaceValue(complex<double> *target, complex<double> *am, int r, double eL, double eR, int order) {

  if(order == 1)
    *target =  1.0/( 1.0*(eR+eL))*( eR*(am[r-1]) + eL*(am[r+1]) );

  if(order == 2)
    *target = -1.0/( 3.0*(eR+eL))*( eR*(-4.0*am[r-1]  + am[r-2]) + eL*( -4.0*am[r+1] + am[r+2]) );

  if(order == 3)
    *target = +1.0/(11.0*(eR+eL))*( eR*(18.0*am[r-1]  - 9.0*am[r-2] + 2.0*am[r-3]) + eL*(18.0*am[r+1] -9.0*am[r+2] + 2.0*am[r+3]) );
}



void WriteEField(double *cx, complex<double> *am, int nx, int r1, int r2, double ncore, const char *fname) {
  int x;

  FILE *gout;
  gout = fopen(fname,"w");

  complex<double> aux = am[nx/2]/(ncore*ncore);

  double nraux = 1.0;
  for(x=0;x<nx;x++) {
    fprintf(gout,"%g %g\n", cx[x], abs(am[x]/aux)/nraux );
    if( x==r1 ) {
      nraux = ncore*ncore;
      fprintf(gout,"%g %g\n", cx[x], abs(am[x]/aux)/nraux );
    }
    if( x==r2 ) {
      nraux = 1.0;
      fprintf(gout,"%g %g\n", cx[x], abs(am[x]/aux)/nraux );
    }
  }
  fclose(gout);
}



int RHS_function(double z, const double y[], double dy[], void *p) {
  RHSpars *pars = (RHSpars*) p;
  int      nc  = pars->nx;
  double   dr  = pars->dr;
  double   la  = pars->la;
  double   nr  = pars->nr;
  double   k01 = 2.0*M_PI/la;
  double   k02 = k01*nr;
  int      iiL = pars->iiL;
  int      iiR = pars->iiR;

  double kef = k02;


  complex<double>  iprop1  = complex<double>(0.0, k01- kef);
  complex<double>  iprop2  = complex<double>(0.0, k02- kef);
  complex<double>  idelta1 = complex<double>(0.0, 1.0/(2*k01*dr*dr) );
  complex<double>  idelta2 = complex<double>(0.0, 1.0/(2*k02*dr*dr) );
  double           dxcoef  = 1.0/2.0/dr;

  complex<double> *am    = (complex<double>*)  y;
  complex<double> *da    = (complex<double>*) dy;


  int r;  
  for(r=1;    r< iiL;  r++) {
    da[r] = idelta1*(am[r-1] -2.0*am[r] + am[r+1])*AA[r] + BB[r]*(am[r+1] - am[r-1])*dxcoef;
    da[r]+= iprop1*am[r];
  }

  for(r=iiL+1;r < iiR;r++) {
    da[r] = idelta2*(am[r-1] -2.0*am[r] + am[r+1]);
    da[r]+= iprop2*am[r];
  }

  for(r=iiR+1;r<nc-1; r++) {
    da[r] = idelta1*(am[r-1] -2.0*am[r] + am[r+1])*AA[r] + BB[r]*(am[r+1] - am[r-1])*dxcoef;
    da[r]+= iprop1*am[r];
  }


  da[0]    = 0.0;
  da[nc-1] = 0.0;

  int iorder = 3;

  complex<double> DC;
  // left interface location
  r  = iiL;
  SetInterfaceValue(&DC, am, r, 1.0, nr*nr, iorder);
  da[r-1] =  idelta1*(am[r-2] -2.0*am[r-1] + DC);
  da[r-1] += iprop1*am[r-1];

  da[r+1] =  idelta2*(DC      -2.0*am[r+1] + am[r+2]);
  da[r+1] += iprop2*am[r+1];

  // it is an 'unused' variable at interface
  da[r-0] = 0.0;


  // right interface location
  r   = iiR;
  SetInterfaceValue(&DC, am, r, nr*nr, 1.0, iorder);
  da[r-1] =  idelta2*(am[r-2] -2.0*am[r-1] + DC);
  da[r-1] += iprop2*am[r-1];

  da[r+1] =  idelta1*(DC      -2.0*am[r+1] + am[r+2]);
  da[r+1] += iprop1*am[r+1];

  // it is an 'unused' variable at interface
  da[r+0] = 0.0;

  return(0);
}


complex<double> modeS(double nb, double r, double cr, double la, double nc, double no) 
{
  // to initialize beam close but not exactly in the waveguide mode

  double be = 2.0*M_PI*nb/la;
  double ko = 2.0*M_PI*no/la;
  double kc = 2.0*M_PI*nc/la;

  double U = sqrt(kc*kc - be*be);
  double W = sqrt(be*be - ko*ko);

  if(r==0.0){
    cout <<"initial field params: " << endl;
    cout << U <<" "<< W <<" "<< kc <<" "<< cr <<" "<< U*sin(U*cr) <<"  * "<< nc*nc/(no*no)*cos(U*cr)*W << endl;
  }

  if( r <= cr ) {
    return(nc*nc*cos(U*r));
  }
  else {
    return(nc*nc*cos(U*cr)*exp(-(r-cr)*W) );
  }
}





int main(int argc, char *argv[]) {
  int x;

  Solver  solver;
  RHSpars rpars;

  /* variables/values as in the matlab scritp */
  double    la     = 1.3e-06;          
  double    nc     = sqrt(2.0);
  double    wc     = 0.75e-06;


  double    LX     = 20.0e-06;
  double    z      = 0.0;
  double    dz     = 1.0e-7;
  double    LZ     = 5.0e-3;

  // reference index
  double    nb;
  nb = 1.30425; 


  /* derived parameters */
  double k0   = 2*M_PI/la;
  double dx   = LX/NX;

  /* to send parameters to RHSfunction */
  rpars.nx  = NX;
  rpars.dr  = dx;
  rpars.la  = la;
  rpars.nr  = nc;
  int r1;
  int r2;
  r1 = rpars.iiL = NX/2 - (int) (0.5*wc/dx);
  r2 = rpars.iiR = NX/2 + (int) (0.5*wc/dx);


  /* coordinates */
  double cx[NX];
  for(x=0;x<NX;x++) cx[x] = dx*(x-NX/2.0 - 0.0);
  cout << "core boundary indices and coordinates: " << endl;
  cout << "r1 : " << r1 <<" "<< cx[r1] << endl;
  cout << "r2 : " << r2 <<" "<< cx[r2] << endl;


  /* initial amplitude */
  complex<double> a0[NX];
  for(x=0;x<NX;x++) a0[x] = modeS(nb, abs(cx[x]), abs(cx[r1+1]), la, nc, 1.0);
  SetInterfaceValue(&a0[r1], a0, r1, 1.0, nc*nc, 3);
  SetInterfaceValue(&a0[r2], a0, r2, nc*nc, 1.0, 3);
  WriteEField(cx, a0, NX, r1, r2, nc, "Xini.dat");

  /* working amplitude */
  complex<double> am[NX];
  for(x=0;x<NX;x++) am[x] = a0[x];

  /* intialize ODE solver */
  solver.Init("rkf45",
             1.0e-3,
             1.0e-3,
             &rpars,
             NX*2,
             RHS_function,
             0,
             true );



/************************** initailize PML ******************************************************/
  int n0;
  complex<double> I = complex<double>(0.0,1.0);
  double ww = 1.0/( LX/8.0 );

  /* these values = no PML */
  for(x=0;x<NX;x++) AA[x] = 1.0;
  for(x=0;x<NX;x++) BB[x] = 0.0;


  /* 30 percent to PML */
  n0 = floor(0.70*NX);
  for(x=n0;x<NX;x++) {
    AA[x] = -1.0/pow(I - 4.0*pow(ww*(cx[x] - cx[n0]),3.0),2.0);
    BB[x] = +6.0*I*ww*ww*ww*pow(cx[x]-cx[n0],2.0)/(k0*pow(I - 4.0*pow(ww*(cx[x] - cx[n0]),3.0),3.0));
  }

  /* 30 percent to PML */
  n0 = floor(0.30*NX);
  for(x=0;x<n0;x++) {
    AA[x] = -1.0/pow(I + 4.0*pow(ww*(cx[x] - cx[n0]),3.0),2.0);
    BB[x] = +6.0*I*ww*ww*ww*pow(cx[x]-cx[n0],2.0)/(k0*pow(I + 4.0*pow(ww*(cx[x] - cx[n0]),3.0),3.0));
  }

  /*************************  end PML adition *****************************************************/



  FILE *fout;
  fout = fopen("LogIMap.dat","w");

  /* evolve amplitude */
  int k,n;  
  double currentz = 0.0;
  double chunk    = 100.0;

  for(n=0;n< chunk;n++) {

    z = 0.0;
    for(k=0;z<LZ/chunk;k++) {
      solver.Step_adaptive(&z,LZ/chunk,&dz,(double*) am);
    }

    for(x=0;x<NX;x++) fprintf(fout,"%g ", log10(norm(am[x]) + 1.0e-20 ));
    fprintf(fout,"\n");

    currentz += z;
    cout << currentz <<" "<< dz  << endl;

    SetInterfaceValue(&am[r1], am, r1, 1.0, nc*nc, 3);
    SetInterfaceValue(&am[r2], am, r2, nc*nc, 1.0, 3);
    WriteEField(cx, am, NX, r1, r2, nc, "Xdim.dat");
  }
  fclose(fout);

  WriteEField(cx, am, NX, r1, r2, nc, "Xdim.dat");

  return(0);
}

