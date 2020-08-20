#ifndef __parameter_io_h__
#define __parameter_io_h__

#include<complex>
#include<list>
#include<vector>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<string>
#include<stdlib.h>

using namespace std;

typedef complex<double> comple;

class Named {

 protected:

  string name;  
  Named* parent;

 public:

  Named() { name = "Unnamed";  parent = 0;}
  
  virtual ~Named() {}

  virtual void Read( istream &stream) {
    if( !stream ) {
      cerr << "PROGSTOP: wrong stream while reading Named:" << name << endl;
      exit(1);
    }

    char   buff[128];
    string saux = "empty_string";

    stream  >> saux;
    while( saux.c_str()[0] == '#' ) {
      stream.getline(buff,128),
      stream  >> saux;
    }

    if( saux != name ) {
      cerr <<"PROGSTOP: named parameter input:\n";
      cerr <<"expected: " << name << endl;
      cerr <<"read    : " << saux << endl;
      exit(1);
    }
  }
  
  virtual void Write(ostream &stream) {
    int lev = Level();
    for(int i=0; i<2*lev; i++) stream <<" ";
    stream << setw(25) << left << name <<" ";
  }

  string &GetName() {return(name);}
  void SetName(const char *given_name) { name = given_name;}
  void SetParent(Named* nm) { parent = nm;}
  int  Level() { if(parent==0) return(0); else return(1 + parent->Level());}
};



class NamedParameterList : public Named {

 protected:

 public:

  vector<Named*> lnp;           

  NamedParameterList() {}

  void IncludeParameter(Named &nm, const char *pname=0) { 
    if(pname) nm.SetName(pname);
    nm.SetParent(this); 
    lnp.push_back(&nm); 
  } 
 
  void Read( istream &stream) {
    Named::Read(stream);
    for( vector< Named* >::iterator ni = lnp.begin(); ni != lnp.end(); ni++ ) (*ni)->Read(stream);
  }

  void Write(ostream &stream) {
    stream << endl;
    Named::Write(stream);
    stream << endl;
    
    for( vector< Named* >::iterator ni = lnp.begin(); ni != lnp.end(); ni++ ) (*ni)->Write(stream);
  }

  void Read(const char *filename) {
    ifstream f;
    f.open(filename);
    if(!f) {
      cout <<"PROGSTOP: can`t open file: " << filename << endl;
      exit(1);
    }
    Read(f);
    f.close();
  }


  Named* FindNamed(const char *name) {
    NamedParameterList *npl;

    for( vector< Named* >::iterator ni = lnp.begin(); ni != lnp.end(); ni++ ) {
      if( (*ni)->GetName() == name ) return(*ni);
      npl = 0;
      npl = dynamic_cast<NamedParameterList*>(*ni);
      if( npl ) npl->FindNamed(name);
    }
    return(0);
  }

};


template <class T, int K = 1 > 
class NamedParameter : public Named {

 public:

  T value[K];

  const string& n()      { return(name    ); }
  T&            v()      { return(value[0]); }
  T&            v(int i) { return(value[i]); }

  void Read( istream &stream) { Named::Read( stream); for(int i=0;i<K;i++) stream >> value[i];  }
  void Write(ostream &stream) { Named::Write(stream); for(int i=0;i<K;i++) stream << value[i] <<" "; stream << endl; }
};


/**************************************************************************/


// the following templates are for bulding parameter hierarchies from inputs: 

// this function should be defined by the framework
template<class T>
T* NameIdentifiedClassPtr_std(const char *name);

// this function should be defined by the user
template<class T>
T* NameIdentifiedClassPtr_usr(const char *name);

// use this macro to build the above functions
#define _getVirtualClassPtr(s,NAME,T) if((string) s == NAME) return( new T ) 

// this is a vector holder for virtual lists
template <class T, int C = 2>
class NamParListVec : public NamedParameterList {

 private:
  
  int length;
  T  *item[C];

 public:
  
  NamedParameter<int>    nc;      // current number of items
  NamedParameter<string> cid[C];  // class ids


  NamParListVec() {
    SetName("FewItems");
    IncludeParameter(nc,"num._of_items");
    length   = 0; 
    for(int i=0;i<C;i++) item[i]=0;
 }

  ~NamParListVec() {
    for(int i=0;i<C;i++) if(item[i]!=0) delete item[i];
  }

  T* Item(int n) {
    if((n<C)&&(item[n]!=0)) return(item[n]); else return(0);
  } 

  int NumOfItems() { return(length); }

  void Read( istream &stream) { 

    NamedParameterList::Read(stream); 

    if(nc.v() > C) {  
      cout <<"Too many items in the list:" << endl;
      //      NamedParameterList::Write(stream); 
      cout <<"Abort" << endl;
      exit(1);
    }

    // read elements
    for(int i=0;i<nc.v();i++) {
      cid[i].SetParent(this);
      cid[i].SetName("type_id_string");
      cid[i].Read(stream);

      if(item[i]==0) item[i] = NameIdentifiedClassPtr_std<T>(cid[i].v().c_str());    
      if(item[i]==0) item[i] = NameIdentifiedClassPtr_usr<T>(cid[i].v().c_str());    

      if(item[i] == 0 ) {
	cout <<"Failed to read or allocate this:" << endl;
	cid[i].Write(cout);
	cout <<"Abort" << endl;
	exit(1);
      }

      item[i]->SetParent(&cid[i]);
      item[i]->Read( stream ); 
      length++;
    }
  }
  
 void Write(ostream &stream) { 
   NamedParameterList::Write(stream); 
   stream << endl;
   for(int i=0;i<nc.v();i++) {
     cid[i].Write(stream);
     item[i]->Write(stream); 
     stream << endl;
   }
  }
};


/**************************************************************************/


// let us give names to freqently used types:
typedef NamedParameter<string>  namstr;
typedef NamedParameter<double>  namdbl;
typedef NamedParameter<int>     namint;

#endif // __parameter_io_h__


