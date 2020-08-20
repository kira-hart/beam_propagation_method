/* this program illustrates implementation and a rudimentary test of DHT */
/* this c-version is more visual w.r.t. what Bessel zeros go where...    */

#include<iostream>
#include<gsl/gsl_sf_bessel.h>

using namespace std;

class DHT {

private:

  int     N;
  double  Rmax;
  double *alphas;
  double *matrix;
  double *coordk;
  double *coordr;

public:

  DHT() {
   alphas = 0;
   matrix = 0;
   coordk = 0;
   coordr = 0;
   N      = 0;
 }

 ~DHT() {
    if(alphas!=0) delete[] alphas;
    if(coordr!=0) delete[] coordr;
    if(coordk!=0) delete[] coordk;
    if(matrix!=0) delete[] matrix;
  }

  int Init(int n, double rmax) {
    alphas = new double[n+1];
    coordk = new double[n];
    coordr = new double[n];
    matrix = new double[n*n];

    if( (alphas==0)||(matrix==0)||(coordk==0)||(coordr==0) ) return(0);

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

	// this part of the matrix is obviously symmetric w.r.t. i<->j
	matrix[j + i*N]  = gsl_sf_bessel_J0(alphas[j]*alphas[i]/alphas[N]);

	// index of the zero that appears in J1 is the second index of the matrix!
        matrix[j + i*N] *= 2.0/(alphas[N]*gsl_sf_bessel_J1(alphas[j])*gsl_sf_bessel_J1(alphas[j]));
      }
    return(N);
  }

  const double *cr() {return(coordr);}
  const double *ck() {return(coordk);}
  const double *TM() {return(matrix);}

  void Transform(double *src, double *tar) {
    for(int i=0;i<N;i++) {
      tar[i] = 0;
      for(int j=0;j<N;j++) tar[i] += matrix[i*N +j]*src[j];
    }
  }

};

int main(int argc, char *argv[]) {
  int i;

  if(argc!=3) {
    cout <<"usage: t.out N m" << endl;
    return(1);
  }

  // this is the order of the Hankel transform
  int M = atoi(argv[1]);

  // Hankel transform object
  DHT dht;
  dht.Init(M,1.0);

  double *vi = new double[M];  // initial vector
  double *sp = new double[M];  // its spatial spectrum
  double *vf = new double[M];  // final vector

  // input vector representing spectrum corresponding to a Bessel mode
  // only one spatial frequency has non-zero amplitude
  for(i=0;i<M;i++) vi[i] = 0.0;

  int mode = atoi(argv[2]);   // this input selects the mode to be 'excited'
  vi[mode] = 1.0;

  // uncomment this for the opposite transformation direction
  // this input should result in a delta-function spectrum
  // for(i=0;i<M;i++) vi[i] = gsl_sf_bessel_J0( dht.ck()[mode]*dht.cr()[i] );

  // two successive transformation should produce the original array
  dht.Transform(vi,sp);
  dht.Transform(sp,vf);

  // show the result: array index, initial vector minus final vector,
  // initial vector, and its spectrum
  // the second column should contain near-zeros
  for(i=0;i<M;i++) {
    cout << i <<" "<< vi[i] -vf[i] <<" "<< sp[i] <<" "<< vf[i] << endl;
  }

  return(0);
}
