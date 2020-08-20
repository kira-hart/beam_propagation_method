/*
Crank-Nicolson BPM illustrated - radially symmetric solution case
Illustrates: Kerr-induced self-focusing collapse
*/

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<complex>

using namespace std;
typedef complex<double> comple;



// for initial condition
comple gaussr(double wr, double r, double z, double k, double f);

// tri-diagonal matrix solver
void TDMsolver(int n, const comple *a, comple *b, const comple *c, comple *v, comple *x);

// to write amplitude array
void Save(const char *fname, int n, double dr, comple *a);

// to measure beams zero and second moments  
void ShapeReport(const comple *a, int n, double dr, double &moment0, double &moment2);


class CNPropagator {

private:
  // dimension
  int     N;

  // LP matrix, RHS, and auxiliary holders
  comple *A;
  comple *B;
  comple *C;
  comple *R;
  comple *X;

  // electric field amplitude holders
  comple *E1;
  comple *E2;

  // nonlinear response 
  comple *N1;
  comple *N2;
  comple  nlcoef;

  // step auxiliaries
  comple  idelta;
  comple *currentE;    // points where current solution resides (either E1 or E2)
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
    X = new comple[N*9];
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
    N1 = &X[N*7];  // nonlinear response holder - current step
    N2 = &X[N*8];  // nonlinear response holder - previous step

    currentE   = E1;
    auxiliaryE = E2;
  }

  void SetIDelta(double dz, double dr, double lambda, double n2, double Iunit) {
    int row;

    // set idelta coeff
    double k0 = 2.0*M_PI/lambda;
    idelta = comple(0.0,  dz/(4.0*k0*dr*dr) );

    // set nonlinear coupling coefficient
    nlcoef = comple(0.0, k0*n2*Iunit*dz);
    for(row=1;row<N;row++) N2[row] = 0.0;

    // construct DELTA operator first
    // the first row correponds to the center r=0, different from the rest
    A[0] =  0.0; // not used, irrelevant
    B[0] = -4.0;
    C[0] = +4.0;
    
    // include radial "correction" to discrete Laplacian:
    for(row=1;row<N;row++) {
      A[row] = +1.0 - 0.5/((double) row);   
      B[row] = -2.0;
      C[row] = +1.0 + 0.5/((double) row);
    }
    
    // define LP matrix using DELTA: LP = 1 -idelta*DELTA
    for(row=0;row<N;row++) {
      A[row] = -idelta*A[row];
      B[row] = +1.0 - idelta*B[row];
      C[row] = -idelta*C[row];
      X[row] = B[row];
    }
  }


  void ExecuteStep() {
    int row;

    // 1. calculate RHS: 

    // diagonal contribution
    for(row=0; row<N; row++) R[row] = currentE[row];

    // first row is special
    R[0] += 4.0*idelta*(currentE[1] - currentE[0]);

    // middle rows
    for(row=1;row<N-1; row++) R[row] += idelta*( currentE[row-1] -2.0*currentE[row] + currentE[row+1] + 0.5*(currentE[row+1] - currentE[row-1])/((double) row) );

    // last row is also special
    row = N-1;
    R[row] += idelta*( currentE[row-1] -2.0*currentE[row] + 0.0 + 0.5*(0.0 - currentE[row-1])/((double) row) );


    /********/

    // RHS modifcation due to nonlinearity
    for(row=0;row<N; row++) N1[row] = nlcoef*currentE[row]*norm(currentE[row]);

    // extrapolation of the nonlinear response -> central point of C-N stencil
    for(row=0;row<N; row++) R[row] = R[row] +  1.5*N1[row] - 0.5*N2[row];

    // save the nonlinear response profile
    for(row=0; row<N; row++) N2[row] = N1[row];

    /********/


    // 2. apply LP inverse

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
  int r;

  // computational domain params
  double R  = 2.0e-02;
  double nr = 5000;
  double dr = R/nr;

  // initial condition params
  double la = 800.0e-09;
  double w0 = 0.5e-02;
  double f  = 0.0;

  //double Iunit = 1.0e+16;
  double Iunit = 1.0e+16;
  double I0    = 1.0;

  // nonlinear index
  double n2 = 1.0e-23;

  // driver params
  double zfinal   = 15.0e-00;
  double dz       = 1.0e-04;
  double zcurrent = 0.0;

  // my solver
  CNPropagator cnp;

  // allocate memory
  cnp.Allocate(nr);

  // set step (and with it the idelta CN coeff)
  cnp.SetIDelta(dz, dr, la, n2, Iunit);

  // set initial condition
  for(r=0;r<nr;r++) cnp.GetAmplitudePtr()[r] = I0*gaussr(w0, dr*r, 0.0, 2.0*M_PI/la, f);
  Save("sol-i.dat",nr,dr,cnp.GetAmplitudePtr());
  
  // evolution
  double zreport = 0.0;
  for(;zcurrent<zfinal;) {
    cnp.ExecuteStep();
    zcurrent += dz;

    if( zcurrent >= zreport) {
      double moment0, moment2;
      ShapeReport(cnp.GetAmplitudePtr(),nr,dr,moment0,moment2);
  
      cout << zcurrent <<" "<< norm(cnp.GetAmplitudePtr()[0]) <<" "<< moment0 <<" "<< moment2/moment0 << endl;
      zreport = zcurrent + 0.001;
    }
    
    if(  norm(cnp.GetAmplitudePtr()[0]) > 3000.0 ) break;
  }

  // final numerical solution
  Save("sol-f.dat",nr,dr,cnp.GetAmplitudePtr());

  // analytic solution to compare with
  for(r=0;r<nr;r++) cnp.GetAmplitudePtr()[r] = gaussr(w0, dr*r, zcurrent, 2.0*M_PI/la, 0.0);
  Save("sol-a.dat",nr,dr,cnp.GetAmplitudePtr());
}




complex<double> gaussr(double wr, double r, double z, double k, double f) 
{
  complex<double> aux1, aux2, icko(0.0,1.0);

  if( wr <= 0.0 ) return(1.0);

  if(f != 0.0) 
    aux1 =  1.0/(wr*wr) + icko*k/(2.0*f);
  else
    aux1 =  1.0/(wr*wr);

  aux2 =  1.0/(1.0 + 2.0*icko*z/k*aux1);
  return( aux2*exp(-r*r*aux1*aux2) );
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


void Save(const char *fname, int n, double dr, comple *a) {
  ofstream f;
  f.open(fname);
  for(int r=0;r<n;r++) f << dr*r <<" "<< norm(a[r]) <<" "<< real(a[r]) <<" "<< imag(a[r]) << endl;
  f.close();
}


void XShapeReport(const comple *a, int n, double dr, double &moment0, double &moment2) {
  int i;
  double r1, r2;

  // this is a rather rough estimate of moments: one should use proper integration weights
  moment0 = 0.0;
  moment2 = 0.0;
  for(i=0;i<n;i++) {
    double r = dr*i;
    moment0 += 2.0*M_PI*dr*norm(a[i])*r;
    moment2 += 2.0*M_PI*dr*norm(a[i])*r*r*r;
  }
}

// this function needs checking!!!
void ShapeReport(const comple *a, int n, double dr,  double &moment0, double &moment2) {
  int i;
  double r1, r2;

  for(i=0;i<n;i++) {
    if(i==0) {
      r2 = dr*0.5;
      
      moment0 = M_PI*norm(a[0])*pow(r2 , 2.0)*0.5;
      moment2 = M_PI*norm(a[0])*pow(r2 , 4.0)*0.25;
    }
    else {
      r1 = dr*(i-0.5);
      r2 = dr*(i+0.5);
      
      moment0 += M_PI*norm(a[i])*(r2*r2 - r1*r1);
      moment2 += M_PI*norm(a[i])*(r2*r2 - r1*r1)*pow((r1+r2)*0.5,2.0);
    }
  }
}
