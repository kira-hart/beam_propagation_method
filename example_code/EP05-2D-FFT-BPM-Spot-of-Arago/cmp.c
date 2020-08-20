/* 
   This is written as to achieve as close corespondence 
   to the Matlab script pSpotOfArago.m as possible.
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<complex.h>
#include<fftw3.h>

/* initial condition = wide gaussian - narower slit */
fftw_complex IC(double r) {

double w0 = 2.0e-03;
double w1 = 5.0e-03;

 if( fabs(r) < w0 ) return(0.0);
 return( exp(-pow(r*r/(w1*w1),4.0) ) );
} 


int main(int argc, char *argv[]) {
  /* array indices */
  int x,y;

  /* for timing    */
  long c0,c1;

  /* variables/values as in the matlab scritp */
  const int NX     = 4096;
  double    lambda = 633.0e-09;          
  double    LX     = 2.0e-02;
  double    dz     = 0.05; 
  double    LZ     = 1.0;

  /* derived parameters */
  double k0   = 2*M_PI/lambda;
  double dk   = 2*M_PI/LX;
  double dx   = LX/NX;
  int    stps = LZ/dz;

  /* allocate amplitude array */
  fftw_complex* am;
  am = fftw_malloc(NX*NX*sizeof(fftw_complex));

  /* allocate spatial spectrum array */
  fftw_complex* sp;
  sp = fftw_malloc(NX*NX*sizeof(fftw_complex));

  /* prepare for fft2d */
  c0 = clock();
  printf("preparing fftw plans\n");

  fftw_plan fftwpF;
  fftw_plan fftwpB;  
  fftwpF = fftw_plan_dft_2d(NX,NX,am,sp,FFTW_FORWARD, FFTW_MEASURE);
  fftwpB = fftw_plan_dft_2d(NX,NX,sp,am,FFTW_BACKWARD,FFTW_MEASURE);

  c1 = clock();
  printf ("\tCPU time spent in plan preparation:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);


  /* initialize arrays */
  c0 = clock();

  /* coordinates */
  double cx[NX];
  for(x=0;x<NX;x++) cx[x] = dx*(x-NX/2);

  /* transverse wavenumbers */
  double kx[NX];
  for(x=0;x<NX/2;x++) kx[x] = dk*x;
  for(   ;x<NX;  x++) kx[x] = dk*(x-NX);

  /* 2D propagator holder */
  fftw_complex* pxy;
  pxy = fftw_malloc(NX*NX*sizeof(fftw_complex));
  for(x=0;x<NX;x++) {
    for(y=0;y<NX;y++) {
      /* paraxial propagator */
      pxy[x+NX*y] = cexp(-I*(kx[x]*kx[x] + kx[y]*kx[y])/(2*k0)*dz )/(NX*NX);
      }
  }
  
  /* 2D boundary guard holder */
  double* bxy;
  bxy = fftw_malloc(NX*NX*sizeof(double));
  for(x=0;x<NX;x++) {
    for(y=0;y<NX;y++) {
      bxy[x+NX*y] = exp(-pow((cx[x]*cx[x] + cx[y]*cx[y])/(LX*LX/5.0),8.0) );
      }
  }
  
  /* amplitude holder, define an initial condition */
  for(x=0;x<NX;x++) {
    for(y=0;y<NX;y++) {
      am[x+NX*y] = IC(sqrt(cx[x]*cx[x] + cx[y]*cx[y]));
    }
  }

  c1 = clock();
  printf ("\tCPU time spent in array initialization:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);


  /* execute steps */
  c0 = clock();
  for(y=1;y<=stps;y++) {
    printf("%d of %d\n",y,stps);
    fftw_execute(fftwpF);
    for(x=0;x<NX*NX;x++) sp[x] *= pxy[x];
    fftw_execute(fftwpB);
    for(x=0;x<NX*NX;x++) am[x] *= bxy[x];
  }

  c1 = clock();
  printf ("\tCPU time spent in steps:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);


  /* save results */
  FILE *fout;
  fout = fopen("amslice.dat","w");
  for(x=0;x<NX;x++) {
    fprintf(fout,"%g %g\n",cx[x],cabs(am[x+NX*(NX/2-1)]));
  }
  fclose(fout);

  /* clean up */
  fftw_destroy_plan(fftwpF);
  fftw_destroy_plan(fftwpB);

  /* should clean up other arrays ... */

  return(0);
}

