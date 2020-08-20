#include<iostream>
#include<math.h>

using namespace std;


int main(int argc, char *argv[]) 
{
  double n2, lambda, waist, intensity, focus;

  cout <<"what n2: ";
  cin >> n2;
  cout << n2 << endl;

  cout <<"what wavelength: ";
  cin >> lambda;
  cout << lambda << endl;

  cout <<"what beam waist: " ;
  cin >> waist;
  cout << waist << endl;

  cout <<"what intensity: ";
  cin >> intensity;
  cout << intensity << endl;

  cout <<"what focal length: ";
  cin >> focus;
  cout << focus << endl;

   double pi,pc,cm,z0,zc,nn;
   
   pi = M_PI/2.0*intensity*pow(waist,2.0);
   pc = pow(0.61*lambda,2.0)*M_PI/(8.0*n2);
   cm = 2.725*sqrt(pow( sqrt(pi/pc) - 0.852 ,2.0) - 0.0219);
   z0 = M_PI*pow(waist,2.0)/lambda;
   
   if( focus != 0.0 )
     zc = z0/(cm + z0/focus);
   else
     zc = z0/cm;
 
   cout << endl; 
   cout << "  #  SF distance: " << endl;
   cout << endl; 
   cout << "  #  Pi:  " << pi << endl;
   cout << "  #  Pc:  " << pc << endl;
   cout << "  #  z0:  " << z0 << endl;
   cout << "  #  zf1: " << zc << endl;
   
   // estimate according to Liu and Chin, Optics Express 13(2005)5754
   double zsf;
   zsf = 0.367 * z0/ sqrt( pow( sqrt(pi/pc) - 0.852, 2.0) - 0.0219 );
   
   double fnl;
   if( fabs(focus) > 0.0 )
     fnl = focus*zsf/(focus + zsf);
   else 
     fnl = zsf;
   
   cout << "  #  zf2: " << fnl << endl;

  return(0);
}
