
#include <iostream>
#include "residue-and-atom-specs.hh"

std::ostream& coot::operator<< (std::ostream& s, const coot::atom_spec_t &spec) {

   s << "[spec: ";
   s << "model ";
   s << spec.model_number;
   s << " ";
   s << "\"";
   s << spec.chain_id;
   s << "\" ";
   s << spec.res_no;
   s << " ";
   s << "\"";
   s << spec.ins_code;
   s << "\"";
   s << " ";
   s << "\"";
   s  << spec.atom_name;
   s << "\"";
   s << " ";
   s << "\"";
   s << spec.alt_conf;
   s << "\"]";

   return s;

}

std::ostream& coot::operator<< (std::ostream& s, const coot::residue_spec_t &spec) {

   if (!spec.unset_p()) { 

      s << "[spec: ";
      // s << "{{debug:: mmdb::MinInt4 is " << MinInt4 << "}} ";
      if (spec.model_number == mmdb::MinInt4)
	 s << "mmdb::MinInt4";
      else
	 s << spec.model_number;
      
      s << " \"";
      s << spec.chain_id;
      s << "\" ";
      s << spec.res_no;
      s << " ";
      s << "\"";
      s << spec.ins_code;
      s << "\"]";
   } else {
      s << "{residue-spec-not-set}";
   } 
   return s;

}


bool
coot::atom_spec_t::matches_spec(mmdb::Atom *atom) const {

   if (atom_name == std::string(atom->name)) {

      if (alt_conf == std::string(atom->altLoc)) {

	 mmdb::Residue *residue_p = atom->residue;
	 
	 if (residue_p) { 
	    
	    if (res_no == atom->GetSeqNum()) {
	       
	       if (ins_code == std::string(atom->GetInsCode())) { 
		  
		  mmdb::Chain *chain_p= atom->GetChain();
		  if (chain_p) {
		     if (chain_id == chain_p->GetChainID()) {
			// std::cout << atom_name << "a complete match " << std::endl;
			return 1;
		     } else {
			// std::cout << atom_name << "a chain mismatch " << std::endl;
			return 0;
		     }
		  } else {
		     // std::cout << atom_name << "a no chain match " << std::endl;
		     // no chain
		     return 1;
		  }
	       } else {
		  // std::cout << atom_name << "an inscode mismatch " << std::endl;
		  return 0;
	       }
	    } else {
	       // std::cout << atom_name << "a resno mismatch " << std::endl;
	       return 0;
	    }
	    
	 } else {
	    // no residue
	    // std::cout << atom_name << "a no chain match " << std::endl;
	    return 1;
	 }
      } else {
	 // std::cout << atom_name << "an altloc mismatch " << std::endl;
	 return 0;
      } 
   } else {
      // std::cout << atom_name << "an atom name mismatch :" << atom->name << ":" << std::endl;
      return 0;
   }
   std::cout << atom_name << " should not happen (matches_spec()) " << atom->name << ":" << std::endl;
   return 0;
}

// return an atom selection handle for the selection in the mol
// that matches the spec.
//
int
coot::residue_spec_t::select_atoms(mmdb::Manager *mol, int selhnd,
				   mmdb::SELECTION_KEY selection_key) {

   if (mol) { 
      mol->SelectAtoms(selhnd, 0, chain_id.c_str(),
		       res_no, ins_code.c_str(),
		       res_no, ins_code.c_str(),
		       "*", "*", "*", "*", selection_key);
   }
   return selhnd;
} 



// the header for this is (in) residue-and-atom-specs.hh.  Hmm... should be fixed.
//
// model_p is a default argument, default 0/NULL (model_number is not set in the atom specs)
//
std::pair<coot::atom_spec_t, coot::atom_spec_t>
coot::link_atoms(mmdb::Link *link, mmdb::Model *model_p) {

   atom_spec_t a1(link->chainID1, link->seqNum1, link->insCode1, link->atName1, link->aloc1);
   atom_spec_t a2(link->chainID2, link->seqNum2, link->insCode2, link->atName2, link->aloc2);

   if (model_p) {
      int mn = model_p->GetSerNum();
      a1.model_number = mn;
      a2.model_number = mn;
   }

   return std::pair<coot::atom_spec_t, coot::atom_spec_t> (a1, a2);
}
