#include"bpm_icn.h"

void SetPML(comple *AA, comple *BB, double *cx, int nx, bool on);

ICN::ICN() {
  NX = NY = NC = NA = 0;
  data = 0;
  pmls = 0;
  epsr = 0;
  cx   = 0;
}

ICN::~ICN() {
  if(data != 0 ) delete[] data;
  if(pmls != 0 ) delete[] pmls;
  if(epsr != 0 ) delete[] epsr;
  if(cx   != 0 ) delete[] cx;
}


void ICN::Allocate(int nx, int ny, double LX, double LY) {
  NX = nx;
  NY = ny;
  NC = 2;
  NA = 3;
  
  int NN = NX*NY*NC*NA;
  data = new comple[NN];
  pmls = new comple[2*NX + 2*NY];
  epsr = new double[NX*NY];
  cx   = new double[NX + NY];
  
  if((data==0)||(epsr==0)||(cx==0)) {
    cout <<"PROGSTOP: ICN::Allocate failed\n";
    exit(1);
  }
  
  // coordinates
  cy = &cx[NX];
  dx =  LX/((double) NX - 1.0);
  dy =  LY/((double) NY - 1.0);
  for(int x=0;x<NX;x++) cx[x] = -LX/2.0 + dx*x;
  for(int y=0;y<NY;y++) cy[y] = -LY/2.0 + dy*y;
  
  // PMLs
  pml_aa_x = &pmls[0];
  pml_bb_x = &pmls[NX];
  pml_aa_y = &pmls[2*NX + 0];
  pml_bb_y = &pmls[2*NX + NY];
  
  PMLinit(true);
  
  // refractive index
  for(int y=0;y<NY;y++) for(int x=0;x<NX;x++) epsr[index(x,y)] = params.epsrel.EpsilonRelative(cx[x],cy[y]);
  
  // initial condition
  comple *ex = EX(0);
  comple *ey = EY(0);
  
  for(int i=0;i<NN;i++) data[i] = 0.0;
  
  for(int y=1;y<NY-1;y++) for(int x=1;x<NX-1;x++) {
      ex[index(x,y)] = params.incond.InitialCondition(cx[x],cy[y],0);
      ey[index(x,y)] = params.incond.InitialCondition(cx[x],cy[y],1);
    }
}

void ICN::Normalize() {
  int i;
  double aux;

  aux = 0.0;
  for(i=0;i<NX*NY*NC;i++) if(aux < abs(data[i])) aux = abs(data[i]);
  aux = 1.0/aux;
  for(i=0;i<NX*NY*NC;i++)  data[i] *= aux;
}


void ICN::PMLinit(bool on) {
  SetPML(pml_aa_x, pml_bb_x, cx, NX, on);
  SetPML(pml_aa_y, pml_bb_y, cy, NY, on);
}


void SetPML(comple *AA, comple *BB, double *cx, int nx, bool on) {
  complex<double> I = complex<double>(0.0,1.0);
  double ww = 1.0/((cx[nx-1] - cx[0])/8.0 );
  int x,n0;

  for(x=0;x<nx;x++) AA[x] = 1.0;
  for(x=0;x<nx;x++) BB[x] = 0.0;
  
  if(on==false) return;

  /* 15 percent to PML */
  n0 = floor(0.85*nx);
  for(x=n0;x<nx;x++) {
    AA[x] = -1.0/pow(I - 4.0*pow(ww*(cx[x] - cx[n0]),3.0),2.0);
    BB[x] = -12.0*ww*pow(cx[x]-cx[n0],2.0)/(pow(I - 4.0*pow(ww*(cx[x] - cx[n0]),3.0),3.0));
  }
  
  /* 15 percent to PML */
  n0 = floor(0.15*nx);
  for(x=0;x<n0;x++) {
    AA[x] = -1.0/pow(I + 4.0*pow(ww*(cx[x] - cx[n0]),3.0),2.0);
    BB[x] = +12.0*ww*pow(cx[x]-cx[n0],2.0)/(pow(I + 4.0*pow(ww*(cx[x] - cx[n0]),3.0),3.0));
  }
}

