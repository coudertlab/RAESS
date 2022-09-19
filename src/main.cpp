#include <local/forcefield.hpp> // define forcefield parameters
#include <local/sphere.hpp>     // Define sampling methodology
#include <local/boxsetting.hpp> // set rectangular shaped box using periodic boundary conditions

#include <gemmi/cif.hpp>        // file -> cif::Document
#include <gemmi/smcif.hpp>      // cif::Document -> SmallStructure
#include <gemmi/symmetry.hpp>   // Space Group manipulation
#include <gemmi/unitcell.hpp>
#include <unordered_map>
// Set key constants
#define R 8.31446261815324e-3 // kJ/mol/K
#define sqrt_2 1.414213562373095
#define min_factor 1.122462048309373  // 2^(1/6)
#define N_A 6.02214076e23    // part/mol

using namespace std;
namespace cif = gemmi::cif;

int main(int argc, char* argv[]) {
  // Set up Input Variables
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  if(argc == 2 && strcmp(argv[1], "--help")==0){
    printf("[USAGE] %s STRUCT FFIELD TEMPER CUTOFF NSTEPS ADSORB REJECT SRADIUS  \n", argv[0]);
    printf("Parameters need to be put in this exact order (a parsing tool will be added soon)\n");
    printf("    STRUCT = Path to material cif file (e.g. KAXQIL_clean.cif)\n");
    printf("    FFIELD = Path to Raspa type forcefield file (e.g. forcefield.def)\n");
    printf("    TEMPER = Temperature in K used in the Boltzmann average (e.g. 298.0)\n");
    printf("    CUTOFF = Cut-off distance in angstrom for the LJ potential (e.g. 12.0)\n");
    printf("    NSTEPS = Number of sampling points around a sampled atom (e.g. 2000)\n");
    printf("    ADSORB = Atomic element tag since only works for monoatomic adsorbent now (e.g. Xe)\n");
    printf("    REJECT = Rejection condition: relative diameter of the rejection threshold distance (e.g. 0.85)\n");
    printf("    SSSIZE = Sampling Sphere Size parameter is the relative distance\n"
           "             compared to sigma as an atomic radius (e.g. 1.6)\n");
    printf("[TEST] %s structure/KAXQIL_clean.cif forcefield/Dreiding_uff.def 298 16 2000 Xe 0.85 1.6  \n", argv[0]);
    exit(0);
  }

  auto structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  int num_steps = stoi(argv[5]);
  string element_guest_str = argv[6];
  double access_coeff_sq = 0;
  if (argv[7]) {access_coeff_sq = stod(argv[7])*stod(argv[7]);}
  double radius_factor = min_factor;
  if (argv[8]) {radius_factor = stod(argv[8]);}

  // Error catch
  if ( num_steps < 0 || temperature < 0 ) {throw invalid_argument( "Received negative value for the Number of Steps or the Temperature or the Accessibility Coefficient" );}
  if ( access_coeff_sq > 1 ) {throw invalid_argument( "Accessibility Coefficient above 1 (Read the purpose of this coeff)" );}

  // Read Forcefield Infos
  LennardJones::Parameters ff_params;
  if (forcefield_path != "DEFAULT") {
    ff_params.read_lj_from_raspa(forcefield_path);
  }

  // Inialize key variables
  string element_host_str;
  double dist = 0;
  double distance_sq = 0;
  double epsilon = 0;
  double sigma = 0;
  double sigma_sq = 0;
  double sigma_6 = 0;
  double exp_energy = 0;

  // Adsorption infos from forcefield dictionary
  pair<double,double> epsilon_sigma = ff_params.get_epsilon_sigma(element_guest_str, true);
  double epsilon_guest = epsilon_sigma.first;
  double sigma_guest = epsilon_sigma.second;

  // Read Structure cif files
  cif::Document doc = cif::read_file(structure_file);
  cif::Block block = doc.sole_block();
  gemmi::SmallStructure structure = gemmi::make_small_structure_from_block(block);
  // Setup spacegroup using number and reset the images properly (don't use hm notations)
  int spacegroup_number = 1;
  for (const char* tag : {"_space_group_IT_number",
                          "_symmetry_Int_Tables_number"})
    if (const string* val = block.find_value(tag)) {
      spacegroup_number = (int) cif::as_number(*val);
      break;
    }
  const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(spacegroup_number);
  structure.cell.set_cell_images_from_spacegroup(sg);
  vector<gemmi::SmallStructure::Site> unique_sites = structure.sites;
  vector<gemmi::SmallStructure::Site> all_sites = structure.get_all_unit_cell_sites();
  // Cell vector
  double a_x = structure.cell.orth.mat[0][0]; double b_x = structure.cell.orth.mat[0][1]; double c_x = structure.cell.orth.mat[0][2];
  double b_y = structure.cell.orth.mat[1][1]; double c_y = structure.cell.orth.mat[1][2];
  double c_z = structure.cell.orth.mat[2][2];
  // Minimal rectangular box that could interact with atoms within the smaller equivalent rectangluar box
  int n_max = int(abs((cutoff + sigma) / a_x)) + 1; 
  int m_max = int(abs((cutoff + sigma) / b_y)) + 1; 
  int l_max = int(abs((cutoff + sigma) / c_z)) + 1; 

  // Creates a list of sites within the cutoff
  vector<array<double,6>> supracell_sites;
  gemmi::Fractional coord;
  string element_host_str_temp = "X";

  unordered_map<string, int> sym_counts;
  for (auto site: all_sites) {
    element_host_str = site.type_symbol;
    if (element_host_str != element_host_str_temp) {
      epsilon_sigma = ff_params.get_epsilon_sigma(element_host_str, false);
      // Lorentz-Berthelot
      epsilon = sqrt( epsilon_sigma.first * epsilon_guest );
      sigma = 0.5 * ( epsilon_sigma.second+sigma_guest );
    }
    element_host_str_temp = element_host_str;
    ++sym_counts[site.label];
    // neighbor list within rectangular box
    move_rect_box(site.fract,a_x,b_x,c_x,b_y,c_y);
    for (int n = -n_max; (n<n_max+1); ++n){
      for (int m = -m_max; (m<m_max+1); ++m) {
        for (int l = -l_max; (l<l_max+1); ++l) {
          // calculate a distance from centre box
          array<double,6> pos_epsilon_sigma;
          coord.x = site.fract.x + n;
          coord.y = site.fract.y + m;
          coord.z = site.fract.z + l;
          gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(coord));
          pos_epsilon_sigma[0] = pos.x;
          pos_epsilon_sigma[1] = pos.y;
          pos_epsilon_sigma[2] = pos.z;
          pos_epsilon_sigma[3] = epsilon;
          pos_epsilon_sigma[4] = sigma * sigma;
          pos_epsilon_sigma[5] = pos_epsilon_sigma[4] * pos_epsilon_sigma[4] * pos_epsilon_sigma[4];
          supracell_sites.push_back(pos_epsilon_sigma);
        }
      }
    }
  }

  if (sym_counts.size() != unique_sites.size()) {throw invalid_argument( "Can't generate symmetry mapping for unique sites, make sure each atoms have a unique label" );}

  vector<gemmi::Vec3> sphere_distr_vector = generateSphereSpirals(num_steps);

  double mass = 0;
  double boltzmann_energy_lj = 0;

  double sum_exp_energy = 0;
  double inv_distance_6;
  double inv_distance_12;
  double sigma_host;
  double radius;
  double energy_lj;
  vector<array<double,6>> neighbor_sites;
  double area_accessible = 0;

  for ( gemmi::SmallStructure::Site site: unique_sites ) {
    // Get LJ parameters
    element_host_str = site.type_symbol;
    sigma_host = ff_params.get_sigma(element_host_str, false);
    radius = radius_factor * 0.5 * (sigma_guest+sigma_host);
    int sym_count = sym_counts[site.label];

    gemmi::Element el(element_host_str.c_str());
    mass += sym_count * el.weight();
    move_rect_box(site.fract,a_x,b_x,c_x,b_y,c_y);
    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));
    // Cell list pruning to have only the sites that are within (cutoff + radius) of the unique site
    neighbor_sites = {};
    for (array<double,6> pos_epsilon_sigma : supracell_sites) {
      double cutoff_2 = cutoff + radius;
      double delta_x = abs(Vsite.x-pos_epsilon_sigma[0]);
      if (delta_x > cutoff_2) {continue;}
      double delta_y = abs(Vsite.y-pos_epsilon_sigma[1]);
      if (delta_y > cutoff_2) {continue;}
      double delta_z = abs(Vsite.z-pos_epsilon_sigma[2]);
      if (delta_z > cutoff_2) {continue;}
      distance_sq = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
      if (distance_sq > cutoff_2*cutoff_2) {continue;}
      neighbor_sites.push_back(pos_epsilon_sigma);
    }
    // Loop around the sphere surface of the unique site
    int count_acc = 0; // count accessible points
    for (gemmi::Vec3 V: sphere_distr_vector) {
      V *= radius;
      // Lennard Jones interaction energies
      energy_lj = 0;
      gemmi::Vec3 pos_neigh;
      bool free = true;
      for(array<double,6> pos_epsilon_sigma : neighbor_sites) {
        double energy_temp = 0;
        pos_neigh = gemmi::Vec3(pos_epsilon_sigma[0], pos_epsilon_sigma[1], pos_epsilon_sigma[2]);
        distance_sq = (V+Vsite).dist_sq(pos_neigh);
        sigma_sq = pos_epsilon_sigma[4];
        if (distance_sq <= sigma_sq*access_coeff_sq) {
          energy_lj = 1e10;
          free = false;
          break;
        }
	else if (distance_sq < cutoff_sq) {
          epsilon = pos_epsilon_sigma[3];
          sigma_6 = pos_epsilon_sigma[5];
          inv_distance_6 = 1.0 / ( distance_sq * distance_sq * distance_sq );
          inv_distance_12 = inv_distance_6 * inv_distance_6;
          energy_temp = epsilon * sigma_6 * ( sigma_6 * (inv_distance_12 - inv_cutoff_12) - inv_distance_6 + inv_cutoff_6 );
          energy_lj += energy_temp;
        }
        if (energy_temp > 0) { free = false; }
      }
      energy_lj *= 4*R;
      if ( free ) {count_acc++; }
      exp_energy = exp(-energy_lj/(R*temperature)); 
      sum_exp_energy += exp_energy;
      boltzmann_energy_lj += exp_energy*energy_lj;
    }
    area_accessible += radius * radius * sym_count * count_acc;
  }
  double Framework_density = 1e-3 * mass/(N_A*structure.cell.volume*1e-30); // kg/m3
  double enthalpy_surface = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  double henry_surface = 1e-3*sum_exp_energy/(R*temperature)/(unique_sites.size()*num_steps)/Framework_density;    // mol/kg/Pa
  area_accessible *= 1e4 * M_PI / (2 * structure.cell.volume * num_steps); // m2/cm3 // Divided by 8 because of sigma and rmin diff
  string structure_name(structure_file);
  structure_name = structure_name.substr(structure_name.find_last_of("/\\") + 1);
  string::size_type const p(structure_name.find_last_of('.'));
  structure_name = structure_name.substr(0, p);
  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)
  cout << structure_name << "," << enthalpy_surface << "," << henry_surface << "," << area_accessible << "," << elapsed_time_ms*0.001 << endl;
}
