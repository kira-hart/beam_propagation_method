/*
Crank-Nicolson BPM illustrated - one-dimensional linear domain
*/

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<complex>

using namespace std;
typedef complex<double> comple;



// for initial condition
comple gauss1D(double wr, double r, double z, double k, double f);

// tri-diagonal matrix solver
void TDMsolver(int n, const comple *a, comple *b, const comple *c, comple *v, comple *x);

// to write amplitude array
void Save(const char *fname, int n, double dx, comple *a);




class CNPropagator {

private:
  // dimension
  int     N;

  // LP matrix, RHS, and auxiliary holders
  comple *A;         // lower sub diagonal
  comple *B;         // diagonal
  comple *C;         // upper diagonal
  comple *R;         // right-hand-side of the linear system to solve
  comple *X;         // permanent storage for diagonal

  // electric field amplitude holders: we need two of these
  comple *E1;
  comple *E2;

  // step auxiliaries
  comple  idelta;
  comple *currentE;    // pointer to where current solution resides (either E1 or E2)
  comple *auxiliaryE;  // this is where the new solution will be stored (either E1 or E2)

public:

  CNPropagator() {
    A=B=C=R=E1=E2=X=0;
  }

  void ReleaseMemory() {
    if(X !=0) delete[] X;
  }

 ~CNPropagator() {
   ReleaseMemory();
  }
  
  void Allocate(int n) {
    ReleaseMemory();

    N = n;
    X = new comple[N*7];
    if(X==0) {
      cerr <<"CNPropagator: failed to allocate\n";
      exit(1);
    }

    A  = &X[N*1];  
    B  = &X[N*2]; 
    C  = &X[N*3]; 
    R  = &X[N*4];  
    E1 = &X[N*5];  
    E2 = &X[N*6];  

    currentE   = E1;  
    auxiliaryE = E2;  
  }

  void SetIDelta(double dz, double dx, double lambda) {

    // set idelta coeff
    double k0 = 2.0*M_PI/lambda;
    idelta = comple(0.0,  dz/(4.0*k0*dx*dx) );

    // construct DELTA operator first
    int row;
    
    for(row=0;row<N;row++) {
      A[row] = +1.0;
      B[row] = -2.0;
      C[row] = +1.0;
    }
    
    // define LM matrix using DELTA: LM = 1 -idelta*DELTA
    for(row=0;row<N;row++) {
      A[row] =      - idelta*A[row];
      B[row] = +1.0 - idelta*B[row];
      C[row] =      - idelta*C[row];
      X[row] = B[row];
    }
  }


  void ExecuteStep() {
    int row;

    // 1. calculate RHS: 

    // diagonal contribution
    for(row=0; row<N; row++) R[row] = currentE[row];

    // first row is special
    row = 0;
    R[row] += idelta*( 0.0 -2.0*currentE[row] + currentE[row+1] );

    // middle rows
    for(row=1;row<N-1; row++) R[row] += idelta*( currentE[row-1] -2.0*currentE[row] + currentE[row+1] );

    // last row is also special
    row = N-1;
    R[row] += idelta*( currentE[row-1] -2.0*currentE[row] + 0.0 );

    // 2. apply LM inverse

    // renew B (because it is destroyed in TDMsolver)
    for(row=0; row<N; row++) B[row] = X[row];
    TDMsolver(N, A, B, C, R, auxiliaryE);


    // 3. update what is current amplitude holder
    comple *aux = currentE;
    currentE   = auxiliaryE;
    auxiliaryE = aux;
  }

  comple* GetAmplitudePtr() {
    return( currentE );
  }
};







int main(int argc, char *argv[]) {
  int x;

  // computational domain params
  double R  = 2.0e-03;
  double nx = 2000;
  double dx = R/nx;

  // initial condition params
  double la = 800.0e-09;
  double w0 = 100.0e-06;

  // driver params
  double zfinal   = 1.0e-01;
  double dz       = 1.0e-04;
  double zcurrent = 0.0;

  // my solver
  CNPropagator cnp;

  // allocate memory
  cnp.Allocate(nx);

  // set step (and with it the idelta CN coeff)
  cnp.SetIDelta(dz, dx, la);

  // set initial condition
  for(x=0;x<nx;x++) cnp.GetAmplitudePtr()[x] = gauss1D(w0, dx*(x - nx/2), 0.0, 2.0*M_PI/la, 0.0);
  Save("sol-i.dat",nx,dx,cnp.GetAmplitudePtr());
  
  // evolution
  for(;zcurrent<zfinal;) {
    cnp.ExecuteStep();
    zcurrent += dz;
  }

  // final numerical solution
  Save("sol-f.dat",nx,dx,cnp.GetAmplitudePtr());

  // analytic solution to compare with
  for(x=0;x<nx;x++) cnp.GetAmplitudePtr()[x] = gauss1D(w0, dx*(x - nx/2), zcurrent, 2.0*M_PI/la, 0.0);
  Save("sol-a.dat",nx,dx,cnp.GetAmplitudePtr());
}




complex<double> gauss1D(double wr, double x, double z, double k, double f) 
{
  complex<double> aux1, aux2, imagi(0.0,1.0);

  if( wr <= 0.0 ) return(1.0);

  if(f != 0.0) 
    aux1 =  1.0/(wr*wr) + imagi*k/(2.0*f);
  else
    aux1 =  1.0/(wr*wr);

  aux2 =  1.0/(1.0 + 2.0*imagi*z/k*aux1);
  return( sqrt(aux2)*exp(-x*x*aux1*aux2) );
}

// Tri-Diagonal Matrix linear system solver
void TDMsolver (int n, const comple *a, comple *b, const comple *c, comple *v, comple *x)
{
        /**
         * n - number of equations
         * a - suB-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
         * b - the main diagonal
         * c - suP-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
         * v - right hand side
         * x - solution
         * NOTE: array b will be DESTROYED!!!
         */
        for (int i = 1; i < n; i++)
        {
                comple m = a[i]/b[i-1];
                b[i] = b[i] - m*c[i-1];
                v[i] = v[i] - m*v[i-1];
        }
 
        x[n-1] = v[n-1]/b[n-1];
 
        for (int i = n - 2; i >= 0; i--)
                x[i]=(v[i]-c[i]*x[i+1])/b[i];
}


void Save(const char *fname, int n, double dx, comple *a) {
  ofstream f;
  f.open(fname);
  for(int x=0;x<n;x++) f << dx*(x - n/2) <<" "<< norm(a[x]) <<" "<< real(a[x]) <<" "<< imag(a[x]) << endl;
  f.close();
}
