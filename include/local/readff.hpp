#include <vector>
#include <map>
#include <string>
#include <fstream>      // file reading

using namespace std;

// TODO Create a class
string ReplaceString(string subject, const string& search, const string& replace) {
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != string::npos) {
    subject.replace(pos, search.length(), replace);
    pos += replace.length();
  }
  return subject;
}

vector<string> SplitString(string s, string sep){
  vector<string> v;
	string temp = "";
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

// Put it in a private / public class system
map<string, vector<string> >  ReadFF(string forcefield_path) {
  map<string, vector<string> > forcefield_dict;
  vector<string> L;
  ifstream MyFile(forcefield_path);
  if (!MyFile) {
      cerr << "Couldn't open forcefield file.\n";
  }
  string myText;
  while (getline (MyFile, myText)) {
    L.push_back(myText);
  }
  vector <string> columns_values = SplitString(L[6].substr(2), ", ");
  vector<string> forcefieldDefInfo(L.begin() + 7, L.end() - 2);

  for (size_t j = 1; j < columns_values.size(); ++j ) {
    forcefield_dict[columns_values[0]].push_back(columns_values[j]);
  }
  for (size_t i = 0; i < forcefieldDefInfo.size(); ++i ) {
    vector<string> split_row_temp = SplitString(ReplaceString(forcefieldDefInfo[i],"\t"," "), " ");
    for (size_t j = 1; j < split_row_temp.size(); ++j ) {
      forcefield_dict[split_row_temp[0]].push_back(split_row_temp[j]);
    }
  }
  return forcefield_dict;
}

vector<string> get_epsilon_sigma(string element, map<string, vector<string> > forcefield_dict) {
  vector<string> epsilon_sigma = {};
  try {
    vector<string> value = forcefield_dict.at(element);
    epsilon_sigma.push_back(value[1]);
    epsilon_sigma.push_back(value[2]);
  }
  catch (const std::out_of_range&) {
    cout << "In get_epsilon_sigma, Key \"" << element.c_str() << "\" not found" << endl;
  }
  return epsilon_sigma;
}