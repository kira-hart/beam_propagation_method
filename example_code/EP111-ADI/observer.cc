#include<observer.h>

void Save2D(string filename, const comple *src,  int nx, int ny, double (*function)(comple)) {
  ofstream f;
  f.open(filename);
  for(int y=0;y<ny;y++) {
    for(int x=0;x<nx;x++) {
      f << function(src[x + y*nx]) <<" ";
    }
    f << endl;
  }
  f.close();
  f.clear();
}


void Save1D(string filename, const comple *src, const double *c, int stride, int n) {
  ofstream f;
  f.open(filename);
  for(int i=0;i<n;i++) {
    f << c[i] <<" "<< norm(src[i*stride]) <<" "<< real(src[i*stride]) <<" "<< imag(src[i*stride]) <<endl;
  }
  f.close();
  f.clear();
}


double fun_real(comple x) { return(real(x)); }
double fun_imag(comple x) { return(imag(x)); }
double fun_norm(comple x) { return(norm(x)); }
double fun_phas(comple x) { return(atan2(real(x),imag(x))); }
