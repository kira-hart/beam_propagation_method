#include"bpm_icn.h"


int main(int argc, char *argv[]) {
  ICN icn;  
  Params_All &params = icn.params;

  params.Read(argv[1]);
  params.Write(cout);

  icn.Allocate(params.domain.NX.v(), params.domain.NX.v(), params.domain.LX.v(), params.domain.LY.v());

  icn.SetStep(comple(0.0,-1.0)*params.driver.dz.v(), params.driver.la.v());


  double current_z = 0.0;
  int s;
  for(s=0;s<100;s++) {

    icn.Normalize();
    icn.Report(0,s);

    //    icn.Report(1,s);
    //    icn.Report(2,s);
    cout <<"Report: " << s <<" "<< current_z << endl;   

    for(int k=0;k<params.driver.rp.v();k++)  {
      icn.Step();
      current_z += params.driver.dz.v();
    }

    if(current_z >= params.driver.zf.v()) break;
  }

  return(0);
}

