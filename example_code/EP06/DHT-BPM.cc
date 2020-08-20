/* bare bone implementation of DHT-based spectral beam propagator */
/* program set up to produce a test data for spot of Arago        */

#include<iostream>
#include<fstream>
#include<complex>
#include<gsl/gsl_sf_bessel.h>

using namespace std;
typedef complex<double> comple;

class DHT_BPM {

private:
public:

  int     N;         // number of grid poits
  double  Rmax;      // domain dimension in real space
  double *alphas;    // this will hold BesselJ0 zeros
  double *matrix;    // transformation matrix holder
  double *coordk;    // coordinate array in k-space
  double *coordr;    // coordinate arra in real, r-space
  double *bguard;    // boundary guard

  comple *propag;    // propagator array
  comple *am0;       // amplitude holders and auxiliaries
  comple *am1;
  comple *am2;


  int Init(int n, double rmax) {
    alphas = new double[n+1];    // this array needs one more slot to hold the "max zero"
    coordk = new double[n];
    coordr = new double[n];
    bguard = new double[n];
    matrix = new double[n*n];
    propag = new comple[n*4];   

    if( (alphas==0)||(matrix==0)||(coordk==0)||(coordr==0)||(bguard==0)||(propag==0) ) return(0);

    am0 = &propag[1*n];   
    am1 = &propag[2*n];
    am2 = &propag[3*n];

    // save the transform dimensions
    N    = n;
    Rmax = rmax;

    int i,j;
    // initialize holder with Bessel J0 zeros
    for(i=0;i<=N;i++) alphas[i] = gsl_sf_bessel_zero_J0(i+1);

    // calculate coordinates in real and sepctral spaces:
    for(i=0;i<N;i++) coordr[i] = Rmax*alphas[i]/alphas[N];
    for(i=0;i<N;i++) coordk[i] =      alphas[i]/Rmax;

    // this the transformation matrix
    for(i=0;i<N;i++) for(j=0;j<N;j++) {

	// M[j + i*n] is M[i][j] 
	matrix[j + i*N]  = gsl_sf_bessel_J0(alphas[j]*alphas[i]/alphas[N]);

	// second matrix index, j, is the one in J1 
        matrix[j + i*N] *= 2.0/(alphas[N]*gsl_sf_bessel_J1(alphas[j])*gsl_sf_bessel_J1(alphas[j]));
      }
    return(N);
  }

  // handles to use from outside
  const double *cr() {return(coordr);}
  const double *ck() {return(coordk);}
  const double *TM() {return(matrix);}

  // transform from source src into target tar:
  void Transform(comple *src, comple *tar) {
    for(int i=0;i<N;i++) {
      tar[i] = 0.0;

      // for non-c-users: i*N + j stands for M[i][j] with j being the column index
      for(int j=0; j<N; j++) tar[i] += matrix[i*N +j]*src[j];
    }
  }

  void SetStep(double dz, comple k0) {
    // non-paraxial propagator
    for(int i=0;i<N;i++) propag[i] = exp(comple(0.0,1.0)*dz*sqrt(k0*k0 - coordk[i]*coordk[i]));

    // bundary guard (shoul normally incorporate dz, I keep it aligned with the other programs)
    for(int i=0;i<N;i++) bguard[i] = exp(-pow( coordr[i]*coordr[i]/(Rmax*Rmax/1.25), 8.0) );
  }

  void Step() {
    Transform(am0,am1);                        // Forward transformation
    for(int i=0;i<N;i++) am1[i] *= propag[i];  // apply propagator
    Transform(am1,am0);                        // return to real space
    for(int i=0;i<N;i++) am0[i] *= bguard[i];  // apply boundary guard
  }


  DHT_BPM() {
   alphas = 0;
   matrix = 0;
   N      = 0;
 }

 ~DHT_BPM() {
    if(alphas!=0) delete[] alphas;
    if(matrix!=0) delete[] matrix;
    // ... should delete the rest of arrays
  }
};


/* initial condition = wide gaussian - narower slit */
comple IC(double r) {

double w0 = 2.0e-03;
double w1 = 5.0e-03;

 if( fabs(r) < w0 ) return(0.0);
 return( exp(-pow(r*r/(w1*w1),4.0) ) );
} 


void Save(const char *fname, const double *c, const comple *a, int n) {
  ofstream f;
  f.open(fname);
  for(int i=0;i<n;i++) {
    f << c[i] <<" "<< sqrt(norm(a[i])) <<" "<< real(a[i]) <<" "<< imag(a[i]) << endl;
  }
  f.close();
}


int main(int argc, char *argv[]) {
  int i;
  const int N = 2048;

  // paramters are hardcode here to match the other simulation for comparison

  DHT_BPM bp;
  bp.Init(N, 1.5e-02);
  bp.SetStep(0.05, 2.0*M_PI/0.633e-06);

  for(int i=0;i<N;i++) bp.am0[i] = IC( bp.cr()[i] );

  Save("Xini.dat",bp.cr(),bp.am0,   bp.N);   // initail condition file
  Save("Xpro.dat",bp.ck(),bp.propag,bp.N);   // file with the propagator to check

  for(int i=0;i<20;i++) {
    bp.Step();
    cout << i << endl;
  }

  Save("Xfin.dat",bp.cr(),bp.am0,bp.N);      // final propgated field

  return(0);
}
