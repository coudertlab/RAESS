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
#define MAX_EXP 1e3 // exp(-1000) = 1.34*10^434, for argument above we skip calculation

double energy_lj_opt(double &epsilon, double sigma_6, double &inv_distance_6, double &inv_cutoff_6, double &inv_distance_12, double &inv_cutoff_12) {
  return epsilon*sigma_6*( sigma_6 * (inv_distance_12 - inv_cutoff_12) - inv_distance_6 + inv_cutoff_6 );
}

int main(int argc, char* argv[]) {
  // Set up Input Variables
  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
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
    printf("[TEST] %s structure/KAXQIL_clean_14.cif forcefield/UFF.def 298 12 2000 Xe 0.85 1.6  \n", argv[0]);
    printf("[RESULT] KAXQIL_clean_14,-44.6853,0.0262976,510.872,0.0826761\n");
    exit(0);
  }

  auto structure_file = argv[1];
  std::string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  int num_steps = stoi(argv[5]);
  std::string element_guest_str = argv[6];
  double access_coeff_sq = 0;
  if (argv[7]) {access_coeff_sq = stod(argv[7])*stod(argv[7]);}
  double radius_factor = min_factor;
  if (argv[8]) {radius_factor = stod(argv[8]);}
  double surface_limitation = 1.5;
  if (argv[9]) {surface_limitation = stod(argv[9]);}
  double MAX_ENERGY = surface_limitation*R*temperature;
  double beta = 1/(R*temperature);

  // Error catch
  if ( num_steps < 0 || temperature < 0 ) {throw invalid_argument( "Received negative value for the Number of Steps or the Temperature or the Accessibility Coefficient" );}
  if ( access_coeff_sq > 1 ) {throw invalid_argument( "Accessibility Coefficient above 1 (Read the purpose of this coeff)" );}

  // Read Forcefield Infos
  ForceField::Parameters ff_params;
  if (forcefield_path != "DEFAULT") {
    ff_params.read_lj_from_raspa(forcefield_path);
  }

  // Adsorption infos from forcefield dictionary
  double sigma_guest = ff_params.get_sigma(element_guest_str, true);
  vector<vector<double>> FF_parameters = ff_params.generate_cross_parameters(element_guest_str);

  // Read Structure cif files
  gemmi::cif::Document doc = gemmi::cif::read_file(structure_file);
  gemmi::cif::Block block = doc.sole_block();
  gemmi::SmallStructure structure = gemmi::make_small_structure_from_block(block);
  // Setup spacegroup using number and reset the images properly (don't use hm notations)
  int spacegroup_number = 1;
  for (const char* tag : {"_space_group_IT_number",
                          "_symmetry_Int_Tables_number"})
    if (const std::string* val = block.find_value(tag)) {
      spacegroup_number = (int) gemmi::cif::as_number(*val);
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
  int l_max = floor( std::abs( (cutoff + sigma_guest) / c_z ) ) + 1;
  int m_max = floor( std::abs( (cutoff + sigma_guest + std::abs(l_max*c_y)) / b_y ) ) + 1;
  int n_max = floor( std::abs( (cutoff + sigma_guest + std::abs(m_max*b_x) + std::abs(l_max*c_x)) / a_x ) ) + 1;
  // center position used to reduce the neighbor list
  gemmi::Position center_pos = gemmi::Position(a_x/2,b_y/2,c_z/2);
  double large_cutoff = cutoff + center_pos.length() + sigma_guest;

  // Creates a list of sites within the cutoff
  vector<array<double,4>> supracell_sites;
  gemmi::Fractional coord;
  std::string element_host_str_temp = "X";

  double molar_mass = 0;
  unordered_map<std::string, int> sym_counts;
  for (auto site: all_sites) {
    std::string element_host_str = site.type_symbol;
    gemmi::Element el_host(element_host_str.c_str());
    int atomic_number = el_host.ordinal();
    molar_mass += el_host.weight();
    ++sym_counts[site.label];
    // neighbor list within rectangular box
    move_rect_box(site.fract,a_x,b_x,c_x,b_y,c_y);
    for (int n = -n_max; (n<n_max+1); ++n) {
      for (int m = -m_max; (m<m_max+1); ++m) {
        for (int l = -l_max; (l<l_max+1); ++l) {
          array<double,4> pos_epsilon_sigma;
          coord.x = site.fract.x + n;
          coord.y = site.fract.y + m;
          coord.z = site.fract.z + l;
          gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(coord));
          double distance_sq = pos.dist_sq(center_pos); 
          if (distance_sq <= large_cutoff*large_cutoff) {
            pos_epsilon_sigma[0] = pos.x;
            pos_epsilon_sigma[1] = pos.y;
            pos_epsilon_sigma[2] = pos.z;
            pos_epsilon_sigma[3] = atomic_number;
            supracell_sites.push_back(pos_epsilon_sigma);
          }
        }
      }
    }  
  }

  if (sym_counts.size() != unique_sites.size()) {throw std::invalid_argument( "Can't generate symmetry mapping for unique sites, make sure each atoms have a unique label" );}

  vector<gemmi::Vec3> sphere_distr_vector = generateSphereSpirals(num_steps);

  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;  
  double area_accessible = 0;

  for ( gemmi::SmallStructure::Site site: unique_sites ) {
    // Get LJ parameters
    std::string element_host_str = site.type_symbol;
    double sigma_host = ff_params.get_sigma(element_host_str, false);
    double radius = radius_factor * 0.5 * (sigma_guest+sigma_host);
    int sym_count = sym_counts[site.label];
    move_rect_box(site.fract,a_x,b_x,c_x,b_y,c_y);
    gemmi::Vec3 Vsite = gemmi::Vec3(structure.cell.orthogonalize(site.fract));
    // Cell list pruning to have only the sites that are within (cutoff + radius) of the unique site
    vector<array<double,4>> neighbor_sites = {};
    double cutoff_2 = cutoff + radius;
    for (array<double,4> pos_epsilon_sigma : supracell_sites) {
      double delta_x = abs(Vsite.x-pos_epsilon_sigma[0]);
      if (delta_x > cutoff_2) {continue;}
      double delta_y = abs(Vsite.y-pos_epsilon_sigma[1]);
      if (delta_y > cutoff_2) {continue;}
      double delta_z = abs(Vsite.z-pos_epsilon_sigma[2]);
      if (delta_z > cutoff_2) {continue;}
      double distance_sq = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
      if (distance_sq > cutoff_2*cutoff_2) {continue;}
      neighbor_sites.push_back(pos_epsilon_sigma);
    }
    // Loop around the sphere surface of the unique site
    int count_acc = 0; // count accessible points
    for (gemmi::Vec3 V: sphere_distr_vector) {
      V *= radius;
      // Lennard Jones interaction energies
      double energy_lj = 0;
      gemmi::Vec3 pos_neigh;
      bool free = true;
      for(array<double,4> pos_epsilon_sigma : neighbor_sites) {
        double energy_temp = 0;
        pos_neigh = gemmi::Vec3(pos_epsilon_sigma[0], pos_epsilon_sigma[1], pos_epsilon_sigma[2]);
        double distance_sq = (V+Vsite).dist_sq(pos_neigh);
        int atomic_number = pos_epsilon_sigma[3];
        double sigma_sq = FF_parameters[atomic_number][2];
        if (distance_sq <= sigma_sq*access_coeff_sq) {
          free = false;
          break;
        }
	      else if (distance_sq <= cutoff_sq) {
          double epsilon = FF_parameters[atomic_number][0];
          double sigma_6 = FF_parameters[atomic_number][3];
          double inv_distance_6 = 1.0 / ( distance_sq * distance_sq * distance_sq );
          double inv_distance_12 = inv_distance_6 * inv_distance_6;
          energy_lj += energy_lj_opt(epsilon, sigma_6, inv_distance_6,inv_cutoff_6, inv_distance_12, inv_cutoff_12);
          if (energy_lj > MAX_ENERGY) { free = false; }
        }
      }
      energy_lj *= 4*R;
      if ( free ) {
        count_acc++; 
        double exp_energy = exp(-beta*energy_lj); 
        sum_exp_energy += exp_energy;
        boltzmann_energy_lj += exp_energy*energy_lj;
      }
    }
    area_accessible += radius * radius * sym_count * count_acc;
  }
  double Framework_density = 1e-3 * molar_mass/(N_A*structure.cell.volume*1e-30); // kg/m3
  double enthalpy_surface = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  double henry_surface = 1e-3*sum_exp_energy*beta/(unique_sites.size()*num_steps)/Framework_density;    // mol/kg/Pa
  area_accessible *= 1e4 * M_PI / (2 * structure.cell.volume * num_steps); // m2/cm3 // Divided by 8 because of sigma and rmin diff
  std::string structure_name(structure_file);
  structure_name = structure_name.substr(structure_name.find_last_of("/\\") + 1);
  std::string::size_type const p(structure_name.find_last_of('.'));
  structure_name = structure_name.substr(0, p);
  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)
  std::cout << structure_name << "," << enthalpy_surface << "," << henry_surface << "," << area_accessible << "," << elapsed_time_ms*0.001 << std::endl;
}
