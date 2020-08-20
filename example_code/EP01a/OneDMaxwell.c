/* 
Bare bones One-D Maxwell: Yee algorithm
*/

#include<stdio.h>
#include<stdlib.h>
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

double InitialElectricField(double x, double x0) {
/* 
in the BPM context, only the electric field is specified, and it is
assumed that the magnetic component is such that the whole beam 
propagates in the desired direction.
*/

  double la       = 25.0;    // how many points per wavelength
  double k0       = 2.0/la;  // wavenumber
  double w0       = 5.0*la;  // pulse duration in cycles

  return( sin(k0*x)*exp(-pow((x - x0)/w0,2.0)) );
}

int main(int argc, char *argv[]) {

  const int N = 4096;
  double EF[N];
  double HF[N];

  double dtoverdx = 0.5;

  int i;
  for(i=0;i<N;i++) {
    // here we try to orchestrate the magnetic field with the electric
    // such that the intial waveform will propagate to the right
    // ... this is almost working ... but not accurately enough ...
    double x0 = 0.5*((double) N);
    double xe = (double) i;

    EF[i] = +InitialElectricField(xe, x0);
    HF[i] = -InitialElectricField(xe, x0);
  }

  Save(EF,HF,N,"E0.dat","H0.dat");  // initial condition

  for(i=0;i<atoi(argv[1]);i++) OneStep(EF,HF,dtoverdx,N);

  Save(EF,HF,N,"E1.dat","H1.dat");  // final fields

  return(0);
}

