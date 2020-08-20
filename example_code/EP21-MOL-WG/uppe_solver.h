#ifndef __solver_h__
#define __solver_h__

/* wrapper for gsl ODE solvers */

#include <iostream>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

using namespace std;


typedef  int (*functiont)(double t, const double y[], double f[], void *p);
typedef  int (*jacobiant)(double t, const double y[], double * dfdy, double dfdt[], void * params);


class Solver {

 private:

 public:

  const 
  gsl_odeiv_step_type *T;
  gsl_odeiv_step      *s;
  gsl_odeiv_control   *c ;
  gsl_odeiv_evolve    *e ;
  gsl_odeiv_system    sys;

  bool adaptive;

  Solver() {
    T = 0;
    s = 0;
    c = 0;
    e = 0;
    adaptive = false;
  }

  ~Solver() {
    if(s) gsl_odeiv_step_free(s);
    if(c) gsl_odeiv_control_free(c);  
    if(e) gsl_odeiv_evolve_free(e);
  }

  void Init(const char *method,  double tol_abs, double tol_rel, 
            void *pars, int n, functiont func,  jacobiant jac=0, 
            bool adapt = false) {
    string tname = method;

    T = 0;
    if( tname == "rk2"   ) T = gsl_odeiv_step_rk2;
    if( tname == "rk4"   ) T = gsl_odeiv_step_rk4;
    if( tname == "rkck"  ) T = gsl_odeiv_step_rkck;
    if( tname == "rkf45" ) T = gsl_odeiv_step_rkf45;
    if( tname == "rk2imp") T = gsl_odeiv_step_rk2imp;
    if( tname == "rk4imp") T = gsl_odeiv_step_rk4imp;
    if( tname == "rk8pd" ) T = gsl_odeiv_step_rk8pd;
    if( tname == "bsimp" ) T = gsl_odeiv_step_bsimp;
    if( tname == "gear1" ) T = gsl_odeiv_step_gear1;
    if( tname == "gear2" ) T = gsl_odeiv_step_gear2;

    if( T==0 ) {
      std::cerr << "WARNING: solver name not reckognized: using default rk4\n";
      T = gsl_odeiv_step_rk4;
    }

    adaptive = adapt;

    s = gsl_odeiv_step_alloc(T, n);
    // TODO: release unused memory (for fixed stp)
    if( 1 ) {
      c = gsl_odeiv_control_y_new(tol_abs, tol_rel);
      e = gsl_odeiv_evolve_alloc(n);
    }

    sys.function  = func;
    sys.jacobian  = jac;
    sys.dimension = n;
    sys.params    = pars;
  }

  void Step_adaptive(double *t, double t1, double *h, double *y) {
    int status;

    if( adaptive == false ) {
      std::cerr <<"WARNING: unitialized ADAPTIVE step attempted\n";
      return;
    }

    status = gsl_odeiv_evolve_apply (e, c, s, &sys, t, t1, h, y);
    if(status != GSL_SUCCESS) {
      std::cerr <<"WARNING: unsuccessful ODE solver step!\n";
    }
  }


  void Step_fixed(double *t, double t1, double *h, double *y) {
    int status;
    //status = gsl_odeiv_step_apply (s, *t, *h, y, e->yerr, 0, e->dydt_out,&sys);
    status = gsl_odeiv_step_apply (s, *t, *h, y, e->yerr, 0, 0,&sys);
    if(status != GSL_SUCCESS) {
      std::cerr <<"WARNING: unsuccessful ODE solver step!\n";
    }
    *t += *h;
  }

  void Reset() {
    gsl_odeiv_step_reset(s);
    gsl_odeiv_evolve_reset(e);
  }
};


#endif // __solver_h__
