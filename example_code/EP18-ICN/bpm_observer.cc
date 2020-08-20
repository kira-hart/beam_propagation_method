#include"bpm_icn.h"

void ICN::SaveMatrix(const char *fname, int part, comple *src) {
  int x,y;
  ofstream f;

  f.open(fname);
  for(y=0;y<NY;y++) {
    if(part==1) 
      for(x=0;x<NX;x++) f << real(src[index(x,y)]) <<" ";
    else 
      for(x=0;x<NX;x++) f << imag(src[index(x,y)]) <<" ";
    f << endl;
  }
  f.close();
}  


string ReportFileName(string basename, string reportid, int component, int count, int rank) {
  ostringstream ostr;
  ostr << basename <<"_"<<  reportid  <<"_"<< component <<"_"<< count; 
  if( rank > -1 )
    ostr <<"_r_"<< rank ;
  ostr<< ends;

  return(ostr.str());
}


void ICN::Report(int c, int n) {
  string fname;

  comple *ex = EX(c);
  comple *ey = EY(c);

  fname = ReportFileName(params.driver.base.v().c_str(),"EXre",c,n,-1);
  SaveMatrix(fname.c_str(),1,ex);
 
  fname = ReportFileName(params.driver.base.v().c_str(),"EXim",c,n,-1);
  SaveMatrix(fname.c_str(),2,ex);

  fname = ReportFileName(params.driver.base.v().c_str(),"EYre",c,n,-1);
  SaveMatrix(fname.c_str(),1,ey);
 
  fname = ReportFileName(params.driver.base.v().c_str(),"EYim",c,n,-1);
  SaveMatrix(fname.c_str(),2,ey);

  if(n==0) {
    ofstream f;
    fname = ReportFileName(params.driver.base.v().c_str(),"epsr",0,n,-1);
    f.open(fname.c_str());

    for(int y=0;y<NY;y++) {
      for(int x=0;x<NX;x++) f << epsr[index(x,y)] <<" ";
      f << endl;
    }

    f.close();
  }
}


