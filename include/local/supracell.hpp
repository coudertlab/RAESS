#include <gemmi/smcif.hpp> 
using namespace std;

void set_supracell(int* n_max, int* m_max, int* l_max, gemmi::SmallStructure structure, double cutoff, double sigma_guest) {

  // Cell parameters
  double a = structure.cell.a; double b = structure.cell.b; double c = structure.cell.c;
  double deg_rad = M_PI/180;
  double alpha = deg_rad*structure.cell.alpha; double beta = deg_rad*structure.cell.beta; double gamma = deg_rad*structure.cell.gamma;
  // Definition of the enlarged cutoff
  double cell_diag_halved = 0.5 * sqrt( a*a + b*b + c*c + 2*(b*c*cos(alpha) + c*a*cos(beta) + a*b*cos(gamma)) );
  double large_cutoff = cutoff + sigma_guest + cell_diag_halved;
  // Cell vector definition
  double a_x = structure.cell.orth.mat[0][0]; double b_x = structure.cell.orth.mat[0][1]; double c_x = structure.cell.orth.mat[0][2];
  double a_y = structure.cell.orth.mat[1][0]; double b_y = structure.cell.orth.mat[1][1]; double c_y = structure.cell.orth.mat[1][2];
  double a_z = structure.cell.orth.mat[2][0]; double b_z = structure.cell.orth.mat[2][1]; double c_z = structure.cell.orth.mat[2][2];
  // Minimal box setting for a triclinic cell to contain a sphere of radius `large_cutoff`
  *n_max = abs(int(large_cutoff * sqrt((b_y*c_x-b_x*c_y)*(b_y*c_x-b_x*c_y)
                      + (b_z*c_x-b_x*c_z)*(b_z*c_x-b_x*c_z)
                      + (b_z*c_y-b_y*c_z)*(b_z*c_y-b_y*c_z))
                      / (a_z*b_y*c_x - a_y*b_z*c_x - a_z*b_x*c_y
                      + a_x*b_z*c_y + a_y*b_x*c_z - a_x*b_y*c_z))) + 1;
  *m_max = abs(int(large_cutoff * sqrt((c_y*a_x-c_x*a_y)*(c_y*a_x-c_x*a_y)
                      + (c_z*a_x-c_x*a_z)*(c_z*a_x-c_x*a_z)
                      + (c_z*a_y-c_y*a_z)*(c_z*a_y-c_y*a_z))
                      / (b_z*c_y*a_x - b_y*c_z*a_x - b_z*c_x*a_y
                      + b_x*c_z*a_y + b_y*c_x*a_z - b_x*c_y*a_z))) + 1;
  *l_max = abs(int(large_cutoff * sqrt((a_y*b_x-a_x*b_y)*(a_y*b_x-a_x*b_y)
                      + (a_z*b_x-a_x*b_z)*(a_z*b_x-a_x*b_z)
                      + (a_z*b_y-a_y*b_z)*(a_z*b_y-a_y*b_z))
                      / (c_z*a_y*b_x - c_y*a_z*b_x - c_z*a_x*b_y
                      + c_x*a_z*b_y + c_y*a_x*b_z - c_x*a_y*b_z))) + 1;
}