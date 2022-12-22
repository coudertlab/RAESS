#include <iostream>      // printing
#include <fstream>       // file reading
#include <string>
#include <vector>
#include <utility>       // std::pair
#include <gemmi/elem.hpp>      // Atomic Element manipulation
#include <cmath>
#include <cstring>

// Gas mono-atomiques pour les adsorbats
#ifndef LOCAL_FF_HPP_
#define LOCAL_FF_HPP_

using std::vector;

namespace ForceField {

std::string ReplaceString(std::string subject, const std::string& search, const std::string& replace) {
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != std::string::npos) {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
  return subject;
}

vector<std::string> SplitString(std::string s, std::string sep){
  vector<std::string> v;
	std::string temp = "";
	for(int i=0;i<s.length();++i){
		if(s.substr(i,sep.size())==sep){
      if(temp!=""){
			  v.push_back(temp);
      }
			temp = "";
      i=i+sep.size()-1;
		}
		else{
			temp.push_back(s[i]);
		}
	}
  if(temp!=""){
    v.push_back(temp);
  }
  return v;
}

struct Parameters {
  enum class MixingRule {
    LorentzBerthelot,
    Jorgensen
  };
  MixingRule mixing_rule = MixingRule::LorentzBerthelot;
  enum class CutoffRule {
    Shifted,
    Truncated
  };
  CutoffRule cutoff_rule = CutoffRule::Shifted;
  enum class InteractionType {
    LennardJones,
    Buckingham
  };
  InteractionType interaction_type = InteractionType::LennardJones;

  static const int N_host = 120;
  double host_epsilons[N_host] = {
    /*X*/ NAN,
    /*H*/  22.1417, /*He*/ 28.1803, /*Li*/ 12.5805, /*Be*/ 42.7737, 
    /*B*/  90.5795, /*C*/ 52.8381, /*N*/ 34.7222, /*O*/ 30.1932, 
    /*F*/  25.161, /*Ne*/ 21.1352, /*Na*/ 15.0966, /*Mg*/ 55.8574, 
    /*Al*/ 254.126, /*Si*/ 202.294, /*P*/ 153.482, /*S*/ 137.882, 
    /*Cl*/ 114.231, /*Ar*/ 93.0956, /*K*/ 17.6127, /*Ca*/ 119.766, 
    /*Sc*/ 9.56117, /*Ti*/ 8.55473, /*V*/ 8.05152, /*Cr*/ 7.5483, 
    /*Mn*/ 6.54186, /*Fe*/ 6.54186, /*Co*/ 7.04508, /*Ni*/ 7.5483, 
    /*Cu*/ 2.5161, /*Zn*/ 62.3992, /*Ga*/ 208.836, /*Ge*/ 190.72, 
    /*As*/ 155.495, /*Se*/ 146.437, /*Br*/ 126.308, /*Kr*/ 110.708, 
    /*Rb*/ 20.1288, /*Sr*/ 118.257, /*Y*/ 36.2318, /*Zr*/ 34.7222, 
    /*Nb*/ 29.69, /*Mo*/ 28.1803, /*Tc*/ 24.1545, /*Ru*/ 28.1803, 
    /*Rh*/ 26.6706, /*Pd*/ 24.1545, /*Ag*/ 18.1159, /*Cd*/ 114.734, 
    /*In*/ 301.429, /*Sn*/ 285.326, /*Sb*/ 225.946, /*Te*/ 200.281, 
    /*I*/ 170.591, /*Xe*/ 167.069, /*Cs*/ 22.6449, /*Ba*/ 183.172, 
    /*La*/ 8.55473, /*Ce*/ 6.54186, /*Pr*/ 5.0322, /*Nd*/ 5.0322, 
    /*Pm*/ 4.52898, /*Sm*/ 4.02576, /*Eu*/ 4.02576, /*Gd*/ 4.52898, 
    /*Tb*/ 3.52254, /*Dy*/ 3.52254, /*Ho*/ 3.52254, /*Er*/ 3.52254, 
    /*Tm*/ 3.01932, /*Yb*/ 114.734, /*Lu*/ 20.632, /*Hf*/ 36.2318, 
    /*Ta*/ 40.7608, /*W*/ 33.7157, /*Re*/ 33.2125, /*Os*/ 18.6191, 
    /*Ir*/ 36.735, /*Pt*/ 40.2576, /*Au*/ 19.6256, /*Hg*/ 193.74, 
    /*Tl*/ 204.383, /*Pb*/ 333.635, /*Bi*/ 260.668, /*Po*/ 163.546, 
    /*At*/ 142.914, /*Rn*/ 124.798, /*Fr*/ 25.161, /*Ra*/ 203.301, 
    /*Ac*/ 16.6063, /*Th*/ 13.0837, /*Pa*/ 11.0708, /*U*/ 11.0708, 
    /*Np*/ 9.56117, /*Pu*/ 8.05152, /*Am*/ 7.04508, /*Cm*/ 6.54186, 
    /*Bk*/ 6.54186, /*Cf*/ 6.54186, /*Es*/ 6.03864, /*Fm*/ 6.03864, 
    /*Md*/ 5.53542, /*No*/ 5.53542,/*Lr*/ NAN, /*Rf*/ NAN, 
    /*Db*/ NAN, /*Sg*/ NAN, /*Bh*/ NAN, /*Hs*/ NAN,
    /*Mt*/ NAN, /*Ds*/ NAN, /*Rg*/ NAN, /*Cn*/ NAN,
    /*Nh*/ NAN, /*Fl*/ NAN, /*Mc*/ NAN, /*Lv*/ NAN, 
    /*Ts*/ NAN, /*Og*/ NAN,
    /*END*/ 0.0
  };

  double host_sigmas[N_host] = {
    /*X*/ NAN,
    /*H*/ 2.57113, /*He*/ 2.1043, /*Li*/ 2.18359, /*Be*/ 2.44552, 
    /*B*/ 3.63754, /*C*/ 3.43085, /*N*/ 3.26069, /*O*/ 3.11815, 
    /*F*/ 2.99698, /*Ne*/ 2.88918, /*Na*/ 2.65755, /*Mg*/ 2.69141, 
    /*Al*/ 4.00815, /*Si*/ 3.82641, /*P*/ 3.69456, /*S*/ 3.59478, 
    /*Cl*/ 3.51638, /*Ar*/ 3.446, /*K*/ 3.39611, /*Ca*/ 3.02816, 
    /*Sc*/ 2.93551, /*Ti*/ 2.8286, /*V*/ 2.80099, /*Cr*/ 2.69319, 
    /*Mn*/ 2.63795, /*Fe*/ 2.5943, /*Co*/ 2.55866, /*Ni*/ 2.52481, 
    /*Cu*/ 3.11369, /*Zn*/ 2.46155, /*Ga*/ 3.90481, /*Ge*/ 3.81305, 
    /*As*/ 3.7685, /*Se*/ 3.74623, /*Br*/ 3.73197, /*Kr*/ 3.68921, 
    /*Rb*/ 3.66516, /*Sr*/ 3.24376, /*Y*/ 2.98006, /*Zr*/ 2.78317, 
    /*Nb*/ 2.81969, /*Mo*/ 2.71902, /*Tc*/ 2.67091, /*Ru*/ 2.63973, 
    /*Rh*/ 2.60944, /*Pd*/ 2.58272, /*Ag*/ 2.80455, /*Cd*/ 2.53728, 
    /*In*/ 3.97608, /*Sn*/ 3.91283, /*Sb*/ 3.93777, /*Te*/ 3.98232, 
    /*I*/ 4.00904, /*Xe*/ 3.92352, /*Cs*/ 4.02419, /*Ba*/ 3.299, 
    /*La*/ 3.13775, /*Ce*/ 3.16804, /*Pr*/ 3.21258, /*Nd*/ 3.18496, 
    /*Pm*/ 3.16002, /*Sm*/ 3.13596, /*Eu*/ 3.11191, /*Gd*/ 3.00055, 
    /*Tb*/ 3.07449, /*Dy*/ 3.054, /*Ho*/ 3.03707, /*Er*/ 3.02104, 
    /*Tm*/ 3.00589, /*Yb*/ 2.98897, /*Lu*/ 3.24287, /*Hf*/ 2.79831, 
    /*Ta*/ 2.82415, /*W*/ 2.73417, /*Re*/ 2.63171, /*Os*/ 2.7796, 
    /*Ir*/ 2.53015, /*Pt*/ 2.45354, /*Au*/ 2.93373, /*Hg*/ 2.40988, 
    /*Tl*/ 204.383, /*Pb*/ 3.82819, /*Bi*/ 3.89323, /*Po*/ 4.19524, 
    /*At*/ 4.23177, /*Rn*/ 4.24513, /*Fr*/ 4.3654, /*Ra*/ 3.27583, 
    /*Ac*/ 3.09855, /*Th*/ 3.02549, /*Pa*/ 3.05044, /*U*/ 3.0246, 
    /*Np*/ 3.05044, /*Pu*/ 3.05044, /*Am*/ 3.01213, /*Cm*/ 2.96313, 
    /*Bk*/ 2.97471, /*Cf*/ 2.95155, /*Es*/ 2.93907, /*Fm*/ 2.92749, 
    /*Md*/ 2.9168, /*No*/ 2.89364, /*Lr*/ NAN, /*Rf*/ NAN, 
    /*Db*/ NAN, /*Sg*/ NAN, /*Bh*/ NAN, /*Hs*/ NAN,
    /*Mt*/ NAN, /*Ds*/ NAN, /*Rg*/ NAN, /*Cn*/ NAN,
    /*Nh*/ NAN, /*Fl*/ NAN, /*Mc*/ NAN, /*Lv*/ NAN, 
    /*Ts*/ NAN, /*Og*/ NAN,
    /*END*/ 0.0
  };

  // The code only handles monoatomic gas adsorption
  int guest_noble_index[N_host] = {
    /*X*/  0,
    /*H*/  0, /*He*/ 1,
    /*Li*/ 0, /*Be*/ 0, /*B*/  0, /*C*/  0, /*N*/  0,
    /*O*/  0, /*F*/  0, /*Ne*/ 2,
    /*Na*/ 0, /*Mg*/ 0, /*Al*/ 0, /*Si*/ 0, /*P*/  0,
    /*S*/  0, /*Cl*/ 0, /*Ar*/ 3,
    /*K*/  0, /*Ca*/ 0, /*Sc*/ 0, /*Ti*/ 0, /*V*/  0,
    /*Cr*/ 0, /*Mn*/ 0, /*Fe*/ 0, /*Co*/ 0, /*Ni*/ 0,
    /*Cu*/ 0, /*Zn*/ 0, /*Ga*/ 0, /*Ge*/ 0, /*As*/ 0,
    /*Se*/ 0, /*Br*/ 0, /*Kr*/ 4,
    /*Rb*/ 0, /*Sr*/ 0, /*Y*/  0, /*Zr*/ 0, /*Nb*/ 0,
    /*Mo*/ 0, /*Tc*/ 0, /*Ru*/ 0, /*Rh*/ 0, /*Pd*/ 0,
    /*Ag*/ 0, /*Cd*/ 0, /*In*/ 0, /*Sn*/ 0, /*Sb*/ 0,
    /*Te*/ 0, /*I*/  0, /*Xe*/ 5,
    /*Cs*/ 0, /*Ba*/ 0, /*La*/ 0, /*Ce*/ 0, /*Pr*/ 0,
    /*Nd*/ 0, /*Pm*/ 0, /*Sm*/ 0, /*Eu*/ 0, /*Gd*/ 0,
    /*Tb*/ 0, /*Dy*/ 0, /*Ho*/ 0, /*Er*/ 0, /*Tm*/ 0,
    /*Yb*/ 0, /*Lu*/ 0, /*Hf*/ 0, /*Ta*/ 0, /*W*/  0,
    /*Re*/ 0, /*Os*/ 0, /*Ir*/ 0, /*Pt*/ 0, /*Au*/ 0,
    /*Hg*/ 0, /*Tl*/ 0, /*Pb*/ 0, /*Bi*/ 0, /*Po*/ 0,
    /*At*/ 0, /*Rn*/ 6,
    /*Fr*/ 0, /*Ra*/ 0, /*Ac*/ 0, /*Th*/ 0, /*Pa*/ 0,
    /*U*/  0, /*Np*/ 0, /*Pu*/ 0, /*Am*/ 0, /*Cm*/ 0,
    /*Bk*/ 0, /*Cf*/ 0, /*Es*/ 0, /*Fm*/ 0, /*Md*/ 0,
    /*No*/ 0, /*Lr*/ 0, /*Rf*/ 0, /*Db*/ 0, /*Sg*/ 0,
    /*Bh*/ 0, /*Hs*/ 0, /*Mt*/ 0, /*Ds*/ 0, /*Rg*/ 0,
    /*Cn*/ 0, /*Nh*/ 0, /*Fl*/ 0, /*Mc*/ 0, /*Lv*/ 0,
    /*Ts*/ 0, /*Og*/ 0, /*END*/ 0
  };

  static const int N_guest = 7; //modify it accordingly
  double guest_epsilons[N_guest] = {
    /*X*/ NAN,
    /*He*/ 10.9, /*Ne*/ NAN, /*Ar*/ 119.8, /*Kr*/ 166.4, 
    /*Xe*/ 221.0, /*Rn*/ NAN
  };
  double guest_sigmas[N_guest] = {
    /*X*/ NAN,
    /*He*/ 2.64, /*Ne*/ NAN, /*Ar*/ 3.34, /*Kr*/ 3.636, 
    /*Xe*/ 4.1,/*Rn*/ NAN
  };

  inline void read_lj_from_raspa(std::string forcefield_path) {
    vector<std::string> L;
    std::ifstream MyFile(forcefield_path);
    if (!MyFile) {
        std::cerr << "Couldn't open forcefield file.\n";
    }
    std::string myText;
    while (getline (MyFile, myText)) {
      L.push_back(myText);
    }
    vector<std::string> ForcefieldInfo(L.begin() + 7, L.end() - 2);

    for (std::string row : ForcefieldInfo) {
      vector<std::string> split_row_temp = SplitString(ReplaceString(row,"\t"," "), " ");
      std::string element_str = split_row_temp[0];
      char endch = element_str.back();
      if (endch == '_') {
        element_str.pop_back();
        gemmi::Element el(element_str);
        if (strcmp(el.name(),"X")!=0) {
          host_epsilons[el.ordinal()] = stod(split_row_temp[2]);
          host_sigmas[el.ordinal()] = stod(split_row_temp[3]);
        }
      }
      else {
        gemmi::Element el(element_str);
        if (strcmp(el.name(),"X")!=0) {
          guest_epsilons[guest_noble_index[el.ordinal()]] = stod(split_row_temp[2]);
          guest_sigmas[guest_noble_index[el.ordinal()]] = stod(split_row_temp[3]);
        }
      }
    }
  }
  
  inline void print_default_lj_params() {
    typedef const char elname_t[3];
    static constexpr elname_t names[119] = {
        "X",  "H",  "He", 
        "Li", "Be", "B",  "C",  "N",  "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
        "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
        "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
        "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
        "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
        "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
        "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
      };
    std::cout << "EPSILON UFF" << std::endl;
    for (int i=0; i<119; i++) {
      gemmi::Element host_el(names[i]);
      std::cout << "/*" << names[i] << "*/ " << host_epsilons[host_el.ordinal()] << ", ";
      if (i % 4 == 0) {
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
    std::cout << "SIGMA UFF" << std::endl;
    for (int i=0; i<119; i++) {
      gemmi::Element host_el(names[i]);
      std::cout << "/*" << names[i] << "*/ " << host_sigmas[host_el.ordinal()] << ", ";
      if (i % 4 == 0) {
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  }

  inline double get_epsilon(std::string element_str, bool is_guest) {
    gemmi::Element el(element_str.c_str());
    if (is_guest) {
      return guest_epsilons[guest_noble_index[el.ordinal()]];
    }
    else {
      return host_epsilons[el.ordinal()];
    }
  }

  inline double get_sigma(std::string element_str, bool is_guest) {
    gemmi::Element el(element_str.c_str());
    if (is_guest) {
      return guest_sigmas[guest_noble_index[el.ordinal()]];
    }
    else {
      return host_sigmas[el.ordinal()];
    }
  }

  inline std::pair<double,double> get_epsilon_sigma(std::string element_str, bool is_guest) {
    gemmi::Element el(element_str.c_str());

    if (is_guest) {
      return std::make_pair(guest_epsilons[guest_noble_index[el.ordinal()]],guest_sigmas[guest_noble_index[el.ordinal()]]);
    }
    else {
      return std::make_pair(host_epsilons[el.ordinal()],host_sigmas[el.ordinal()]);
    }
  }

  vector<vector<double>> mixing_LJ_LB(std::string element_host) {
    vector<vector<double>> LJ_parameters_mixed;
    gemmi::Element el_host(element_host.c_str());
    double epsilon_guest = guest_epsilons[guest_noble_index[el_host.ordinal()]];
    double sigma_guest = guest_sigmas[guest_noble_index[el_host.ordinal()]];
    for (size_t i_host=0; i_host<N_host; i_host++){
      double epsilon = sqrt(epsilon_guest*host_epsilons[i_host]);
      double sigma = 0.5 * (sigma_guest+host_sigmas[i_host]);
      double sigma_sq = sigma * sigma;
      double sigma_6 = sigma_sq * sigma_sq * sigma_sq;
      vector<double> epsilon_sigma_sq_6{epsilon, sigma, sigma_sq, sigma_6};
      LJ_parameters_mixed.push_back(epsilon_sigma_sq_6);
    }
    return LJ_parameters_mixed;
  }

  vector<vector<double>> mixing_LJ_Jorgensen(std::string element_host) {
    vector<vector<double>> LJ_parameters_mixed;
    gemmi::Element el_host(element_host.c_str());
    double epsilon_guest = guest_epsilons[guest_noble_index[el_host.ordinal()]];
    double sigma_guest = guest_sigmas[guest_noble_index[el_host.ordinal()]];
    for (size_t i_host=0; i_host<N_host; i_host++){
      double epsilon = sqrt(epsilon_guest*host_epsilons[i_host]);
      double sigma = sqrt(sigma_guest*host_sigmas[i_host]);
      double sigma_sq = sigma * sigma;
      double sigma_6 = sigma_sq * sigma_sq * sigma_sq;
      vector<double> epsilon_sigma_sq_6{epsilon, sigma, sigma_sq, sigma_6};
      LJ_parameters_mixed.push_back(epsilon_sigma_sq_6);
    }
    return LJ_parameters_mixed;
  }

  vector<vector<double>> generate_cross_parameters(std::string element_host) {
    if (interaction_type == InteractionType::LennardJones) {
      if (mixing_rule == MixingRule::LorentzBerthelot){
        return mixing_LJ_LB(element_host);
      }
      else if (mixing_rule == MixingRule::Jorgensen){
        return mixing_LJ_Jorgensen(element_host);
      }
      else {throw std::invalid_argument( "This type of mixing rule is not implemented yet" );}
    }
    else {throw std::invalid_argument( "The type of forcefield given is not implemented yet" );}
  }
};

} // LennardJones
#endif