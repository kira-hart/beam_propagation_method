#include"bpm_icn.h"


void ICN::SetStep(comple delta_z, double lambda) {
  comple dzz   = delta_z;
  dz   = abs(delta_z);
  nref = sqrt(params.epsrel.EpsilonRelative(0.0,0.0));
  k0   = 2.0*M_PI/lambda;
  
  if( imag(dzz) != 0.0 ) PMLinit(false);

  coeff_xx = comple(0.0,1.0)*dzz/(2.0*k0*nref*dx*dx); 
  coeff_yy = comple(0.0,1.0)*dzz/(2.0*k0*nref*dy*dy); 
  coeff_di = comple(0.0,1.0)*dzz*k0/(2.0*nref);
  coeff_xy = comple(0.0,1.0)*dzz/(8.0*k0*nref*dx*dy); 

}


void ICN::Add(int src, int tar, double coeff) {
  int i,N;

  comple *S = EX(src);
  comple *T = EX(tar);

  N = NX*NY*2;
  for(i=0;i<N;i++) T[i] += S[i]*coeff;
}

void ICN::ApplyRHS(int src, int tar) {
  ApplyAXX(src,  tar);
  ApplyAYY(src,  tar);
  ApplyAXY(src,  tar);
  ApplyAYX(src,  tar);
}


void ICN::Step() {
  /* formula (13) in the ICN paper */
  ApplyRHS(0,1);
  Add(1, 0, 1.00);

  ApplyRHS(1,2);
  Add(2, 0, 0.50);

  ApplyRHS(2,1);
  Add(1, 0, 0.25);
}

void ICN::ApplyAXX(int src, int tar) {
  int x,y;

  comple *T = EX(tar);
  comple *S = EX(src);

  for(x=0;x<NX;x++) T[index(x,   0)] = 0.0;
  for(x=0;x<NX;x++) T[index(x,NY-1)] = 0.0;
  for(y=0;y<NY;y++) T[index(0,   y)] = 0.0;
  for(y=0;y<NY;y++) T[index(NX-1,y)] = 0.0;
   

  for(y=1;y<NY-1;y++) {
    for(x=1;x<NX-1;x++) {
      double epsL = epsr[index(x-1,y)];
      double epsC = epsr[index(x,  y)];
      double epsR = epsr[index(x+1,y)];
      
      double auxL = (epsL + epsC)/(2.0*epsC);
      double auxC = (epsL + epsC)/(2.0*epsL) + (epsC + epsR)/(2.0*epsR);
      double auxR = (epsR + epsC)/(2.0*epsC);
      
      T[index(x,y)]  =  coeff_xx*pml_aa_x[x]*(auxR*S[index(x+1,y)] - auxC*S[index(x,y)] + auxL*S[index(x-1,y)]);
      
      T[index(x,y)] +=  coeff_yy*pml_aa_y[y]*(     S[index(x,y+1)] -  2.0*S[index(x,y)] +      S[index(x,y-1)]);
      
      T[index(x,y)] +=  coeff_xx*pml_bb_x[x]*( S[index(x+1,y)] - S[index(x-1,y)] );
						  
      T[index(x,y)] +=  coeff_yy*pml_bb_y[y]*( S[index(x,y+1)] - S[index(x,y-1)] );

      T[index(x,y)] +=  coeff_di*(epsC - nref*nref)*S[index(x,y)];
    }
  }
}

void ICN::ApplyAYY(int src, int tar) {
  int x,y;

  comple *T = EY(tar);
  comple *S = EY(src);

  for(x=0;x<NX;x++) T[index(x,   0)] = 0.0;
  for(x=0;x<NX;x++) T[index(x,NY-1)] = 0.0;
  for(y=0;y<NY;y++) T[index(0,   y)] = 0.0;
  for(y=0;y<NY;y++) T[index(NX-1,y)] = 0.0;

  for(y=1;y<NY-1;y++) {
    for(x=1;x<NX-1;x++) {
      double epsL = epsr[index(x,y-1)];
      double epsC = epsr[index(x,  y)];
      double epsR = epsr[index(x,y+1)];
      
      double auxL = (epsL + epsC)/(2.0*epsC);
      double auxC = (epsL + epsC)/(2.0*epsL) + (epsC + epsR)/(2.0*epsR);
      double auxR = (epsR + epsC)/(2.0*epsC);
      
      T[index(x,y)]  =  coeff_xx*pml_aa_x[x]*(     S[index(x+1,y)] -  2.0*S[index(x,y)] +      S[index(x-1,y)]);
      
      T[index(x,y)] +=  coeff_yy*pml_aa_y[y]*(auxR*S[index(x,y+1)] - auxC*S[index(x,y)] + auxL*S[index(x,y-1)]);
      
      T[index(x,y)] +=  coeff_xx*pml_bb_x[x]*( S[index(x+1,y)] - S[index(x-1,y)] );
						  
      T[index(x,y)] +=  coeff_yy*pml_bb_y[y]*( S[index(x,y+1)] - S[index(x,y-1)] );

      T[index(x,y)] +=  coeff_di*(epsC - nref*nref)*S[index(x,y)];
    }
  }
}

void ICN::ApplyAXY(int src, int tar) {
  int x,y;

  comple *T = EX(tar);
  comple *S = EY(src);

  for(x=0;x<NX;x++) T[index(x,   0)] = 0.0;
  for(x=0;x<NX;x++) T[index(x,NY-1)] = 0.0;
  for(y=0;y<NY;y++) T[index(0,   y)] = 0.0;
  for(y=0;y<NY;y++) T[index(NX-1,y)] = 0.0;

  for(y=1;y<NY-1;y++) {
    for(x=1;x<NX-1;x++) {
      double epsPP = epsr[index(x+1,y+1)];
      double epsPM = epsr[index(x+1,y-1)];
      double epsMP = epsr[index(x-1,y+1)];
      double epsMM = epsr[index(x-1,y-1)];
      double epsP0 = epsr[index(x+1,y  )];
      double epsM0 = epsr[index(x-1,y  )];
      
      T[index(x,y)]  +=  (+(epsPP/epsP0 - 1.0)*S[index(x+1,y+1)]   
                          -(epsMP/epsM0 - 1.0)*S[index(x-1,y+1)] 
	                  -(epsPM/epsP0 - 1.0)*S[index(x+1,y-1)]
	                  +(epsMM/epsM0 - 1.0)*S[index(x-1,y-1)])*coeff_xy;
    }
  }
}

void ICN::ApplyAYX(int src, int tar) {
  int x,y;

  comple *T = EY(tar);
  comple *S = EX(src);

  for(x=0;x<NX;x++) T[index(x,   0)] = 0.0;
  for(x=0;x<NX;x++) T[index(x,NY-1)] = 0.0;
  for(y=0;y<NY;y++) T[index(0,   y)] = 0.0;
  for(y=0;y<NY;y++) T[index(NX-1,y)] = 0.0;

  for(y=1;y<NY-1;y++) {
    for(x=1;x<NX-1;x++) {
      double epsPP = epsr[index(x+1,y+1)];
      double epsPM = epsr[index(x+1,y-1)];
      double epsMP = epsr[index(x-1,y+1)];
      double epsMM = epsr[index(x-1,y-1)];
      double eps0P = epsr[index(x,  y+1)];
      double eps0M = epsr[index(x,  y-1)];
    }
  }
}

