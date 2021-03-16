

#ifndef ATOM_PULL_HH
#define ATOM_PULL_HH

#include <clipper/core/coords.h>
#include "geometry/residue-and-atom-specs.hh"

class atom_pull_info_t {

   bool status;
public:
   clipper::Coord_orth pos;
   coot::atom_spec_t spec;

   atom_pull_info_t() { status = false; }
   atom_pull_info_t(const coot::atom_spec_t &spec_in, const clipper::Coord_orth &pos_in) :
      pos(pos_in), spec(spec_in) {
      status = true;
   }

   void off() { status = false; }

   void on() { status = true; }

   bool get_status() const { return status; }

   // return first false if number on not found
   std::pair<bool, int> find_spec(mmdb::PAtom *atoms, int n_atoms) const {

      int idx = -1;
      bool local_status = false;
      if (status) { 
	 for (int iat=0; iat<n_atoms; iat++) {
	    if (coot::atom_spec_t(atoms[iat]).is_same(spec)) {
	       idx = iat;
	       local_status = true;
	       break;
	    } 
	 }
      } 
      return std::pair<bool, int> (local_status, idx);
   }
};

#endif // ATOM_PULL_HH

