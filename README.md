# RAESS
RAESS (Rapid Adsorption Enthalpy Surface Sampling) is a sampling tool for rapid adsorption enthalpy and Henry constant calculations in nanoporous materials 

## Installation

Check if you have `c++11` compiler installed (may work with other compiler but has been tested and mainly used with this compiler).

Tested for c++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0

Compilation that follows the rules set in the `Makefile`:
```
make all
```

## Example Usage

If you want to run a surface sampling simulation on the structure KAXQIL (CSD code) from CoRE MOF 2019 all-solvent removed with the uff forcefield at 298K with a 12A cutoff for the xenon. 
```
./raess structure/KAXQIL_clean_14.cif forcefield/UFF.def 298 12 2000 Xe 0.85 1.6
```
You should get an output that has values close to this:
```
KAXQIL_clean_14,-44.4186,0.0283404,70.9159,0.0707517
```

The results are printed in a comma separated format: structure, adsorption enthalpy (kJ/mol), Henry coefficient (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)

This binary is supposed to be used in a high-throughput manner to add rows to a csv file.

## Paper

This work is presented in the following article (please cite):
[“Rapid Adsorption Enthalpy Surface Sampling (RAESS) to Characterize Nanoporous Materials”](https://doi.org/10.1039/D2SC05810C), E. Ren and F.-X. Coudert, _Chem. Sci._, **2023**, DOI: [10.1039/D2SC05810C](https://doi.org/10.1039/D2SC05810C)

The current implementation is a bit faster than the one presented in the paper while remaining as accurate. (0.314 seconds computation time in average on the CoREMOF 2019 dataset with a LCD higher than 3.7 angstrom instead of 0.34 seconds in the paper).

## Acknowledgement

This code has been developed during a PhD thesis co-financed by the CEA and Orano under the supervision of François-Xavier Coudert: https://github.com/fxcoudert

This code includes the library developed in Gemmi: 
https://github.com/project-gemmi/gemmi.git

## License

The MIT License (MIT)

Copyright (c) 2022 Emmanuel Ren

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
