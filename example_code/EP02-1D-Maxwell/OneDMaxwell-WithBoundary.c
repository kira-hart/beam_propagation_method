/* 
Bare bones One-D Maxwell: Yee algorithm
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

void OneStep(double *E, double *H, double dtoverdx, int n) {
  int i;

  // since we do the update in-place, save these auxiliaries 
  double aux_L_I = E[1];    // this is E_I at left end
  double aux_L_B = E[0];    // this is E_B at left end
  double aux_R_I = E[n-2];  // this is E_I at right end
  double aux_R_B = E[n-1];  // this is E_B at right end

  // update inside of the domain as before
  for(i=1;i<n-1;i++) E[i] = E[i] + dtoverdx*(H[i+1] - H[i  ]);

  // ABCs:
  E[0  ] = aux_L_I + (dtoverdx - 1.0)/(dtoverdx + 1.0)*(E[1]   - aux_L_B );
  E[n-1] = aux_R_I + (dtoverdx - 1.0)/(dtoverdx + 1.0)*(E[n-2] - aux_R_B );

  // update inside of the domain as before
  for(i=1;i<n  ;i++) H[i] = H[i] + dtoverdx*(E[i  ] - E[i-1]);

  // this grid point is not used - it is outside
  H[0]   = 0.0;
}


void Save(double *E, double *H, int n, const char *Ename, const char *Hname) {
  int i;
  FILE *fout;

  fout = fopen(Ename,"w");
  for(i=0;i<n;i++) fprintf(fout,"%g %g\n",i+0.0,E[i]);
  fclose(fout);

  fout = fopen(Hname,"w");
  for(i=1;i<n;i++) fprintf(fout,"%g %g\n",i-0.5,H[i]);
  fclose(fout);
}

double InitialElectricField(double x, double x0) {
/* 
in the BPM context, only the electric field is specified, and it is
assumed that the magnetic component is such that the whole beam 
propagates inthe desired direction.
*/

  double la       = 25.0;

  la /= 2.0;

  double k0       = 2.0/la;
  double w0       = 5.0*la;

  return( sin(k0*x)*exp(-pow((x - x0)/w0,2.0)) );
}

int main(int argc, char *argv[]) {

  const int N = 4096;
  double EF[N];
  double HF[N];

  double dtoverdx = 0.5;

  int i;
  for(i=0;i<N;i++) {
    double x0 = 0.5*((double) N);
    double xe = (double) i;
    double xh = xe - 0.5;

    xh -= 0.5*dtoverdx;

    EF[i] = +InitialElectricField(xe, x0);

    // this will create two pulse in opposite directions
    HF[i] = -0.0*InitialElectricField(xh, x0);
  }

  Save(EF,HF,N,"E0.dat","H0.dat");

  for(i=0;i<atoi(argv[1]);i++) OneStep(EF,HF,dtoverdx,N);

  Save(EF,HF,N,"E1.dat","H1.dat");

  return(0);
}

