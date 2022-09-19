#include <gemmi/smcif.hpp> 
#include <cmath>
using namespace std;

double modulo_1(double a) {
  if (a<0) {return fmod(a,1)+1;}
  else {return fmod(a,1);}
}

// Move the atom into the rectangular box using PBC
void move_rect_box(gemmi::Fractional &Vsite, double a_x, double b_x, double c_x, double b_y, double c_y) {
  // double x = a_x*Vsite.x + b_x*Vsite.y + c_x*Vsite.z;
  // double y =               b_y*Vsite.y + c_y*Vsite.z;
  // double z =                             c_z*Vsite.z;
  Vsite.z = modulo_1(Vsite.z);              // z between 0 and 1
  if (c_y==0) {
    Vsite.y = modulo_1(Vsite.y);
  }
  else {
    double temp = (c_y*Vsite.z)/b_y;
    Vsite.y = modulo_1(Vsite.y+temp) - temp; // y between 0 and 1
  }
  if (b_x==0 && c_x==0) {
    Vsite.x = modulo_1(Vsite.x);
  }
  else {
    double temp = (b_x*Vsite.y + c_x*Vsite.z)/a_x;
    Vsite.x = modulo_1(Vsite.x+temp) - temp; // x between 0 and 1
  }
}
