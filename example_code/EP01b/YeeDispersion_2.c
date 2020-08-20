/* 
Bare bones One-D Maxwell: Yee algorithm

Measurement of the dispersion relation. 
Based on a three-sample formula explained in the notes.
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

void OneStep(double *E, double *H, double dtoverdx, int n) {
  int i;

  for(i=0;i<n-1;i++) E[i] = E[i] + dtoverdx*(H[i+1] - H[i  ]);
  E[n-1] = E[0];

  for(i=1;i<n  ;i++) H[i] = H[i] + dtoverdx*(E[i  ] - E[i-1]);
  H[0]   = H[n-1];
}


void Save(double *E, double *H, int n, const char *Ename, const char *Hname) {
  int i;
  FILE *fout;

  fout = fopen(Ename,"w");
  for(i=0;i<n;i++) fprintf(fout,"%g %g\n",i+0.0,E[i]);
  fclose(fout);

  fout = fopen(Hname,"w");
  for(i=0;i<n;i++) fprintf(fout,"%g %g\n",i-0.5,H[i]);
  fclose(fout);
}

double InitialElectricField(double x, double k0) {
  return( cos(k0*x) );
}

int main(int argc, char *argv[]) {

  const int N = 4096;
  double EF[N];
  double HF[N];

  double dtoverdx = 0.75;

  int el;
  double k0;

  for(el=N/2-1;el>0;el-=200) {
  k0 = el*2.0*M_PI/((double) N - 1.0);    


  int i;
  for(i=0;i<N;i++) {
    double xe = (double) i;

    EF[i] = +InitialElectricField(xe, k0);
    HF[i] = 0.0;
  }


  double EFold = EF[0];
  int    start;
  int    istart, count;

  start  = 0;
  istart = 0;
  count  = 0;

  EFold = EF[0];
  OneStep(EF,HF,dtoverdx,N);
  OneStep(EF,HF,dtoverdx,N);

  double omega  = 2.0/dtoverdx*acos( sqrt(EF[0] + 3.0*EFold)/(2.0*sqrt(EFold)) );
  double yeeval = 2.0/dtoverdx*asin( dtoverdx*sin(k0/2.0));

  printf("%E %E %E\n",k0,omega,yeeval);
  }


  return(0);
}

