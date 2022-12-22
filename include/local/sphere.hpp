#include <random>
#include <vector>
#include <cmath>
#include <chrono>
#include <gemmi/math.hpp>

double invsqrtQuake( double number ) {
  double y = number;
  double x2 = y * 0.5;
  std::int64_t i = *(std::int64_t *) &y;
  // The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
  i = 0x5fe6eb50c7b537a9 - (i >> 1);
  y = *(double *) &i;
  y = y * (1.5 - (x2 * y * y));   // 1st iteration
  return y;
  }

// Direction independent
vector<gemmi::Vec3> generateSphereNormalRandom(int num_steps) {
  vector<gemmi::Vec3> v;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::normal_distribution<double> normal_distrib (0.0,1.0);
  gemmi::Vec3 coord = gemmi::Vec3(0,0,0);
  double norm_sq = 0;
  for (int i = 0; i < num_steps; i++) {
    coord = gemmi::Vec3(normal_distrib(generator), normal_distrib(generator), normal_distrib(generator));
    norm_sq = coord.length_sq();
    v.push_back(coord*invsqrtQuake(norm_sq));
  }
  return v;
}

// Projection of a cylindre on a sphere (cos_phi being the height) Archimedes Theorem
vector<gemmi::Vec3> generateSphereAngleRandom(int num_steps) {
  vector<gemmi::Vec3> v;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  double cos_theta; double sin_theta;
  double cos_phi; double sin_phi;
  gemmi::Vec3 coord;
  for (int i = 0; i < num_steps; i++) {
    double theta = 2 * M_PI * uniform01(generator);
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    cos_phi = 1 - 2 * uniform01(generator);
    sin_phi = sqrt(1 - cos_phi * cos_phi);
    coord = gemmi::Vec3( sin_phi * cos_theta, sin_phi * sin_theta, cos_phi );
    v.push_back(coord);
  }
  return v;
}

vector<gemmi::Vec3> generateSphereCubeRandom(int num_steps) {
  vector<gemmi::Vec3> v;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> uniform01(-1.0, 1.0);
  double x; double y; double z;
  double norm_sq;
  gemmi::Vec3 coord;
  for (int i = 0; i < num_steps; i++) {
    x = uniform01(generator); y = uniform01(generator); z = uniform01(generator);
    norm_sq = x*x + y*y + z*z;
    while ( norm_sq > 1 ) {
      x = uniform01(generator); y = uniform01(generator); z = uniform01(generator);
      norm_sq = x*x + y*y + z*z;
    }
    coord = gemmi::Vec3(x,y,z);
    v.push_back(coord.normalized());
  }
  return v;
}

// Golden number spirals
vector<gemmi::Vec3> generateSphereSpirals(int num_steps) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  vector<gemmi::Vec3> v;
  gemmi::Vec3 coord;
  double d_theta = M_PI * (3-sqrt(5)) ;
  double theta = 2 * M_PI * uniform01(generator);
  double d_cos_phi = 2.0/num_steps;
  double cos_phi = 1 + d_cos_phi*0.5;
  double sin_theta; double cos_theta;
  double sin_phi;
  for (int i = 0; i < num_steps; i++) {
    theta += d_theta;
    sin_theta = sin(theta); 
    cos_theta = cos(theta);
    cos_phi -= d_cos_phi;
    sin_phi = sqrt(1 - cos_phi * cos_phi);
    coord = gemmi::Vec3( sin_phi * cos_theta, sin_phi * sin_theta, cos_phi );
    v.push_back(coord);
  }
  return v;
}
