/* src/atom-name-bits.cc
 * 
 * Copyright 2008 by the University of Oxford
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

namespace coot { 
   // e.g. 
   // "Mg" -> <"  MG", "MG">
   // "I"  -> <"   I", "IOD">
   // 
   class atom_name_bits_t {
   public:
      atom_name_bits_t() { filled = false; }
      bool filled;
      std::string atom_name;
      std::string element_name;
      std::string res_name;
      atom_name_bits_t(const std::string &type) {
	 filled = false;
	 if (type == "Br") {
	    atom_name = "BR  ";
	    element_name = "BR";
	    res_name = "BR";
	    filled = true;
	 }
	 if (type == "Ca") {
	    atom_name = "CA  ";
	    element_name = "CA";
	    res_name = "CA";
	    filled = true;
	 }
	 if (type == "Na") {
	    atom_name = "NA  ";
	    element_name = "NA";
	    res_name = "NA";
	    filled = true;
	 }
	 if (type == "Cl") {
	    atom_name = "CL  ";
	    element_name = "CL";
	    res_name = "CL";
	    filled = true;
	 }
	 if (type == "I") {
	    atom_name = " I  ";
	    element_name = "I";
	    res_name = "IOD";
	    filled = true;
	 }
	 if (type == "Mg") {
	    atom_name = "MG  ";
	    element_name = "MG";
	    res_name = "MG";
	    filled = true;
	 }
	 if (! filled) {
	    // make up (guess) the residue type and element
	    std::string at_name = util::upcase(type);
        atom_name = at_name;
        res_name = at_name;
        element_name = at_name;
	    if (type.length() > 4)
	       atom_name = at_name.substr(0,4);
	    if (type.length() > 3)
	       res_name = at_name.substr(0,3);
	    if (type.length() > 2)
	       element_name = at_name.substr(0,2);
	    filled = true;
	 }
      }
      void SetAtom(mmdb::Atom *at, mmdb::Residue *res) {
	 if (filled) { 
	    at->SetAtomName(atom_name.c_str());
	    at->SetElementName(element_name.c_str());
	    res->SetResName(res_name.c_str());
	 }
      } 
   };
}
