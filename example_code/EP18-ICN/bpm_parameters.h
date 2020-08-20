#include"parameter_input_output.h"

class Params_Driver : public NamedParameterList {
 public:

  namdbl la;
  namdbl dz;
  namdbl zf;
  namint rp;
  namstr base;

  Params_Driver() {
    SetName("Driver");
    IncludeParameter(base,"basename");
    IncludeParameter(la,"wavelength");
    IncludeParameter(dz,"z_step");
    IncludeParameter(zf,"z_stop");
    IncludeParameter(rp,"report_every");
  }

};

class Params_Domain : public NamedParameterList {
 public:

  namdbl LX;
  namdbl LY;
  namint NX;
  namint NY;

  Params_Domain() {
    SetName("Computational_Domain");
    IncludeParameter(LX,"size_LX");
    IncludeParameter(NX,"points_NX");
    IncludeParameter(LY,"size_LY");
    IncludeParameter(NY,"points_NY");
  }

};


class Params_Incond :  public NamedParameterList {
 public:

  namdbl bsize;
  namdbl angle;
  namdbl amplx;
  namdbl amply;

  double lambda;

  Params_Incond() {
    SetName("Initial_condition");
    IncludeParameter(amplx, "E_x");
    IncludeParameter(amply, "E_y");
    IncludeParameter(bsize, "beam_size");
    IncludeParameter(angle, "angle");
  }



  complex<double> gaussx(double wr, double r, double z, double k, double f, double kx) 
    {
      complex<double> aux1, aux2, icko(0.0,1.0);
      
      if( wr <= 0.0 ) return(1.0);
      
      if(f != 0.0) 
	aux1 =  1.0/(wr*wr) + icko*k/(2.0*f);
      else
	aux1 =  1.0/(wr*wr);
      aux2 =  1.0/(1.0 + 2.0*icko*z/k*aux1);
      
      return( sqrt(aux2)*exp(-r*r*aux1*aux2 + complex<double>(0.0,kx*r)));
    }

  
  comple InitialCondition(double rx, double ry, int c) {
    double lambda = 800.0e-09;
    double k0     = 2*M_PI/lambda;
    double kx     = k0*angle.v()/360.0*2*M_PI;
    double w0     = bsize.v();
    
    if( c == 0 )
      return( amplx.v()*gaussx(w0, rx, 0.0, k0, 0.0, +1.0*kx) *gaussx(w0, ry, 0.0, k0, 0.0, +1.0*kx) );
    else
      return( amply.v()*gaussx(w0, rx, 0.0, k0, 0.0, -1.0*kx) *gaussx(w0, ry, 0.0, k0, 0.0, +1.0*kx) );
  }
  
};

class Params_Epsrel :  public NamedParameterList {
 public:

  namdbl ns;
  namdbl nc;
  namdbl d0;
  namdbl d1;
  namdbl d2;
  namdbl d3;

  Params_Epsrel() {
    SetName("Refractive_index");
    IncludeParameter(ns,"substrate");
    IncludeParameter(nc,"core");
    IncludeParameter(d0,"dimension_0");
    IncludeParameter(d1,"dimension_1");
    IncludeParameter(d2,"dimension_2");
    IncludeParameter(d3,"dimension_3");
  }


  double EpsilonRelative(double rx, double ry) {
    if( ry < d0.v() ) return( ns.v()*ns.v() );
    if( ry < d1.v() ) return( nc.v()*nc.v() );
    if( (ry < d2.v())&&(abs(rx) < d3.v()) ) return( nc.v()*nc.v() );

    return( 1.0 );
  }

};



class Params_All : public NamedParameterList {
 public:

  Params_Driver driver;
  Params_Domain domain;
  Params_Incond incond;
  Params_Epsrel epsrel;

  Params_All() {
    SetName("BPM_simulation_parameters");
    IncludeParameter(driver,"Driver");
    IncludeParameter(domain,"Computational_domain");
    IncludeParameter(incond,"Initial_condition");
    IncludeParameter(epsrel,"Refractive_index");
  }

};




