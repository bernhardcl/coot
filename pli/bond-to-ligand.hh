#ifndef COOT_PLI_BOND_TO_LIGAND_HH
#define COOT_PLI_BOND_TO_LIGAND_HH

#include <string>

// contains a coordinate on the flat ligand (input mol coords) that
// corresponds to an atom on the ligand - when we add read the
// ligand mol file, we'll attach to an atom that it needs a bond to
// the residue (that encapsulates this bond description).
//
class bond_to_ligand_t {
private:
   bool is_set_;
public:
   // sync to c-interface-ligands.hh
   enum { H_BOND_DONOR_MAINCHAIN,
      H_BOND_DONOR_SIDECHAIN,
      H_BOND_ACCEPTOR_MAINCHAIN,
      H_BOND_ACCEPTOR_SIDECHAIN,
      METAL_CONTACT_BOND,
      BOND_COVALENT,
      BOND_OTHER };
   std::string ligand_atom_name;
   double bond_length;
   int bond_type; // acceptor or donor
   bond_to_ligand_t(const std::string &n, double b) : ligand_atom_name(n) {
      bond_length = b;
      bond_type = 1;
      is_set_ = 1;
   }
   bond_to_ligand_t() { is_set_ = 0; bond_type = 0; bond_length = 0.0; }
   bool is_set() const { return is_set_; }
};



#endif // COOT_PLI_BOND_TO_LIGAND_HH
