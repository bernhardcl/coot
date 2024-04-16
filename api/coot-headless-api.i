
%module coot_headless_api

%{
#include "molecules-container.hh"
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"


namespace std {
   %template(vector_string) vector<std::string>;
   %template(pairbf) pair<bool, float>;
   %template(IntVector) vector<int>;
}

namespace std {
  %template(coot_chain_validation)   std::vector<coot::chain_validation_information_t>;
  %template(coot_residue_validation) std::vector<coot::residue_validation_information_t>;
}

%init %{
  // init_coot_as_python_module();
%}

%feature("autodoc", "1"); // add doc string for Intellisense (hopefully)

%include "molecules-container.hh"

