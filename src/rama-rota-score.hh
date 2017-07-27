/* src/rama-rota-score.hh
 * 
 * Copyright 2010 by the University of Oxford
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

   // ------------ molecule probability scoring ------------
   class rama_score_t {
   public:
      rama_score_t() {
	 score = 0.0;
	 score_non_sec_str = 0.0;
	 n_zeros = 0;
      }
      // for all residues
      std::vector<std::pair<residue_spec_t, double> >  scores;
      // for non-Secondary structure residues
      std::vector<std::pair<residue_spec_t, double> >  scores_non_sec_str;
      double score;
      double score_non_sec_str;
      int n_residues() const { return scores.size(); }
      int n_residues_non_sec_str() const { return scores_non_sec_str.size(); }
      int n_zeros;
      std::vector<std::pair<residue_spec_t, int> > region;
   };

   // ==-------------- all molecule rotamer scoring --------------
   class rotamer_score_t {
   public:
      rotamer_score_t() {
	 score = 0.0;
	 n_pass = 0;
      } 
      std::vector<std::pair<residue_spec_t, double> > scores;
      double score;
      int n_pass; // GLY, PRO, ALA
      int n_rotamer_residues() const { return scores.size(); }
      void add (const residue_spec_t &rs, double p) {
	 std::pair<residue_spec_t, double> pair(rs,p);
	 scores.push_back(pair);
      }
   };
   

}
