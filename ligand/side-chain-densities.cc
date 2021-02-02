
#include <fstream>
#include <iomanip>
#include <thread>

#include <boost/math/distributions/skew_normal.hpp>

#include "compat/coot-sysdep.h"
#include "analysis/stats.hh"
#include "utils/coot-utils.hh"
#include "side-chain-densities.hh"
#include "richardson-rotamer.hh"

coot::density_box_t::density_box_t(float *density_box_in,
                                   mmdb::Residue *residue_p_in,
                                   int n_steps_in) {
   init();
   density_box = density_box_in;
   residue_p = residue_p_in;
   n_steps = n_steps_in;
   mean=0;
   var = 0;
}


void
coot::side_chain_densities::get_results_addition_lock() {

   bool unlocked = false;
   while (! results_addition_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::microseconds(1));
      unlocked = false;
   }
}

void
coot::side_chain_densities::release_results_addition_lock() {

   results_addition_lock = false;
}

coot::side_chain_densities::side_chain_densities() {

   int n_steps = 5;
   float grid_box_radius = 5.0; // half the width/length of the box (not diagonal)

   std::string n_steps_str = std::to_string(n_steps);
   std::string grid_box_radius_str = util::float_to_string_using_dec_pl(grid_box_radius, 1);
   std::string pdd = package_data_dir(); // xxx/share/coot
   std::string fn = "useable-grid-points-nstep=" + n_steps_str + ",box_radius=" +
      grid_box_radius_str + "-charybdis.data";
   std::string dir_1 = util::append_dir_dir(pdd,   "data");
   std::string dir_2 = util::append_dir_dir(dir_1, "assign-side-chains");
   std::string pathed_file_name = util::append_dir_file(dir_2, fn);
   std::string side_chain_data_sub_dir = "side-chain-data";

   std::string pathed_scd_dir_name =  util::append_dir_file(dir_2, side_chain_data_sub_dir);

   init(n_steps, grid_box_radius, pathed_file_name);
   set_data_dir(pathed_scd_dir_name);

}

void
coot::density_box_t::self_normalize() {

   int n = 2 * n_steps + 1;
   int nnn = n * n * n;
   double sum = 0;
   int n_grid_points = 0;
   for (int i=0; i<nnn; i++) {
      if (density_box[i] > 0.0) {
         sum += density_box[i];
         n_grid_points++;
      }
   }
   if (n_grid_points > 0) {
      double av = sum/static_cast<double>(n_grid_points);
      double sc = 1.0/av;
      for (int i=0; i<nnn; i++)
         if (density_box[i] > -1000)
            density_box[i] *= sc;
   }
}

bool
coot::side_chain_densities::like_the_others(const std::map<int, std::string> &chain,
                                            const std::vector<std::map<int, std::string> > &other_chains) const {

   bool is_like_the_others = false;
   unsigned int n_chain = chain.size();
   for (std::size_t i=0; i<other_chains.size(); i++) {
      const std::map<int, std::string> &other_chain = other_chains[i];
      unsigned int n = n_chain;
      if (other_chains.size() < n)
         n = other_chains.size();
      unsigned int n_match = 0;
      std::map<int, std::string>::const_iterator it;
      for (it=chain.begin(); it!=chain.end(); it++) {
         const int &key = it->first;
         const std::string &val = it->second;

         std::map<int, std::string>::const_iterator it_other = other_chain.find(key);
         if (it_other != other_chain.end())
            if (val == it_other->second)
               n_match++;
      }

      float frac = static_cast<float>(n_match) / static_cast<float>(n);
      if (frac > 0.7) {
         is_like_the_others = true;
         break;
      }
   }

   return is_like_the_others;
}

std::map<int, std::string>
coot::side_chain_densities::make_sequence_for_chain(mmdb::Chain *chain_p) const {

   std::map<int, std::string> m;
   int n_residues = chain_p->GetNumberOfResidues();
   for (int i=0; i<n_residues; i++) {
      mmdb::Residue *residue_p = chain_p->GetResidue(i);
      if (residue_p) {
         int res_no = residue_p->GetSeqNum();
         std::string res_name = residue_p->GetResName();
         m[res_no] = res_name;
      }
   }
   return m;

}

void coot::side_chain_densities::proc_mol(const std::string &id,
                                          mmdb::Manager *mol,
                                          const clipper::Xmap<float> &xmap) {

   // don't sample chains that have a sequence similar to chains that we've
   // processed before
   //
   std::vector<std::map<int, std::string> > done_chains;

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         if (chain_p) {
            std::map<int, std::string> sequence_for_chain = make_sequence_for_chain(chain_p);
            if (! like_the_others(sequence_for_chain, done_chains))
               proc_chain(id, chain_p, xmap);
            done_chains.push_back(sequence_for_chain);
         }
      }
   }

   normalize_density_boxes(id);
   write_density_boxes();
   for (std::size_t i=0; i<density_boxes.size(); i++)
      density_boxes[i].clear();
}

void
coot::side_chain_densities::fill_useable_grid_points_vector(const std::string &file_name) {

   if (! file_name.empty()) {
      std::ifstream f(file_name.c_str());
      if (f) {
         std::string line;
         while (std::getline(f, line)) {
            std::vector<std::string> words = coot::util::split_string_no_blanks(line);
            if (words.size() == 1) {
               int idx = util::string_to_int(words[0]);
               useable_grid_points.insert(idx);
            }
         }
      } else {
         std::cout << "ERROR:: side_chain_densities::fill_useable_grid_points_vector file name not found "
                   << file_name << std::endl;
      }
   }

}

bool
coot::side_chain_densities::is_close_to_atoms(const std::vector<std::pair<double, clipper::Coord_orth> > &atom_positions,
                                              const clipper::Coord_orth &test_position) const {

   bool is_close = false;
   for (std::size_t i=0; i<atom_positions.size(); i++) {
      double max_atom_radius_sqrd = atom_positions[i].first;
      double d = (atom_positions[i].second - test_position).lengthsq();
      if (d < max_atom_radius_sqrd) {
         is_close = true;
         break;
      }
   }
   return is_close;

}

std::string coot::side_chain_densities::get_rotamer_name(mmdb::Residue *res) const {
   std::string res_name = res->GetResName();
   std::string alt_conf;
   mmdb::Manager *mol = 0;
   richardson_rotamer rr(res, alt_conf, mol, 0.0, 1);
   rotamer_probability_info_t prob = rr.probability_of_this_rotamer();
   std::string rotamer_name = util::remove_whitespace(prob.rotamer_name);
   return rotamer_name;

}

clipper::Coord_orth
coot::side_chain_densities::make_pt_in_grid(int ix, int iy, int iz, const float &step_size,
                                            const std::vector<clipper::Coord_orth> &axes) const {

   clipper::Coord_orth pt(0,0,0);
   pt += ix * step_size * axes[0];
   pt += iy * step_size * axes[1];
   pt += iz * step_size * axes[2];

   return pt;

}

// return the "guessed" sequence
std::string
coot::side_chain_densities::probability_of_each_rotamer_at_each_residue(mmdb::Manager *mol,
                                                                        const std::string &chain_id,
                                                                        int resno_start, int resno_end,
                                                                        const clipper::Xmap<float> &xmap,
                                                                        bool verbose_output_mode) {

   std::vector<std::pair<mmdb::Residue *, std::string> > best_guess; // a bit of fun

   bool all_chain = false;
   if (resno_end < resno_start) all_chain = true;

   // What is the probability of each rotamer at each residue?
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         if (chain_p) {
            std::string this_chain_id(chain_p->GetChainID());
            if (this_chain_id == chain_id) {
               int n_residues = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_residues; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     std::string res_name = residue_p->GetResName();
                     bool do_it = all_chain;
                     if (res_no >= resno_start) {
                        if (res_no <= resno_end) {
                           do_it = true;
                        }
                     }

                     if (do_it) {
                        std::map<std::string, double> likelihood_map =
                           likelihood_of_each_rotamer_at_this_residue(residue_p, xmap, false, verbose_output_mode);
                        std::map<std::string, double>::const_iterator it;
                        double best_score = -999999999999999.9;
                        std::string best_type;
                        for (it=likelihood_map.begin(); it!=likelihood_map.end(); it++) {
                           const std::string &res_name = it->first;
                           const double &score = it->second;
                           if (score > best_score) {
                              best_score = score;
                              best_type = res_name;
                           }
                        }
                        std::pair<mmdb::Residue *, std::string> p(residue_p, best_type);
                        best_guess.push_back(p);
                     }
                  }
               }
            }
         }
      }
   }

   std::string best_guess_sequence;
   for (std::size_t i=0; i<best_guess.size(); i++) {
      mmdb::Residue *res = best_guess[i].first;
      std::string best_guess_for_residue = best_guess[i].second;
      if (verbose_output_mode) // this is for devel/analysis/testing, where the answer is known
         std::cout << res->GetChainID() << " " << res->GetSeqNum()
                   << " real-type: "  << res->GetResName()
                   << " best-guess: " << best_guess_for_residue << std::endl;
      best_guess_sequence += single_letter_code(best_guess_for_residue);
   }

   // std::cout << "best guess sequence:\n" << best_guess_sequence << std::endl;

   return best_guess_sequence;
}

std::map<std::string, double>
coot::side_chain_densities::likelihood_of_each_rotamer_at_this_residue(mmdb::Residue *residue_p,
                                                                       const clipper::Xmap<float> &xmap,
                                                                       bool limit_to_correct_rotamers_only,
                                                                       bool verbose_output_mode) {
   return get_rotamer_likelihoods(residue_p, xmap, limit_to_correct_rotamers_only, verbose_output_mode);
}

char
coot::side_chain_densities::single_letter_code(const std::string &res_name) const {

   char r = '.';
   if (res_name == "ALA") r = 'A';
   if (res_name == "ARG") r = 'R';
   if (res_name == "ASN") r = 'N';
   if (res_name == "ASP") r = 'D';
   if (res_name == "CYS") r = 'C';
   if (res_name == "GLN") r = 'Q';
   if (res_name == "GLY") r = 'G';
   if (res_name == "GLU") r = 'E';
   if (res_name == "HIS") r = 'H';
   if (res_name == "ILE") r = 'I';
   if (res_name == "LEU") r = 'L';
   if (res_name == "LYS") r = 'K';
   if (res_name == "MET") r = 'M';
   if (res_name == "PHE") r = 'F';
   if (res_name == "PRO") r = 'P';
   if (res_name == "SER") r = 'S';
   if (res_name == "THR") r = 'T';
   if (res_name == "TRP") r = 'W';
   if (res_name == "TYR") r = 'Y';
   if (res_name == "VAL") r = 'V';

   return r;
}

std::vector<mmdb::Residue *>
coot::side_chain_densities::make_a_run_of_residues(mmdb::Manager *mol, const std::string &chain_id,
                                                   int resno_start, int resno_end) const {

   std::vector<mmdb::Residue *> a_run_of_residues;

   // What is the probability of each rotamer at each residue?
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         if (chain_p) {
            std::string this_chain_id(chain_p->GetChainID());
            if (this_chain_id == chain_id) {
               int n_residues = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_residues; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no = residue_p->GetSeqNum();
                     if (res_no >= resno_start) {
                        if (res_no <= resno_end) {
                           a_run_of_residues.push_back(residue_p);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return a_run_of_residues;
}
void
coot::side_chain_densities::test_sequence(mmdb::Manager *mol,
                                          const std::string &chain_id, int resno_start, int resno_end,
                                          const clipper::Xmap<float> &xmap,
                                          const std::string &sequence_name, // from fasta file
                                          const std::string &sequence) {
   std::vector<results_t> results;
   std::string gene_name = sequence_name;
   std::vector<std::string> parts = util::split_string_no_blanks(gene_name);
   if (parts.size() > 0) gene_name = parts[0];

   // What is the probability of each rotamer at each residue?

   std::vector<mmdb::Residue *> a_run_of_residues = make_a_run_of_residues(mol, chain_id, resno_start, resno_end);

   std::cout << "debug:: in test_sequence() with " << a_run_of_residues.size()
             << " residues in a_run_of_residues" << std::endl;

   if (! a_run_of_residues.empty()) {
      fill_residue_blocks(a_run_of_residues, xmap); // return fast if already filled.
      int n_residues = a_run_of_residues.size();
      std::vector<std::pair<mmdb::Residue *, std::map<std::string, double> > > scored_residues(n_residues);
      for (int i=0; i<n_residues; i++) {
         mmdb::Residue *residue_p = a_run_of_residues[i];
         std::map<std::string, double> likelihood_map =
            likelihood_of_each_rotamer_at_this_residue(residue_p, xmap);
         std::pair<mmdb::Residue *, std::map<std::string, double> > p(residue_p, likelihood_map);

         if (false) { // debugging
            std::cout << "debug:: transfer to scored_residues " << i<< " " << residue_spec_t(residue_p) << " "
                      << residue_p->GetResName() << " " << " with score map: ---" << std::endl;
            std::map<std::string, double>::const_iterator it_debug;
            for (it_debug=likelihood_map.begin();
                 it_debug!=likelihood_map.end();
                 it_debug++) {
               std::cout << "   " << it_debug->first << " " << it_debug->second
                         << std::endl;
            }
            std::cout << "-----" << std::endl;
         }

         scored_residues[i] = p;
         // std::cout << "debug:: scored_residues - storing map of size " << p.second.size() << " for residue vector index "
         // << i << " " << residue_spec_t(p.first) << std::endl;
      }

      std::string true_sequence;
      for (int i=0; i<n_residues; i++) {
         mmdb::Residue *residue_p = a_run_of_residues[i];
         std::string res_name(residue_p->GetResName());
         char s = single_letter_code(res_name);
         true_sequence += s;
      }

      // slide
      std::cout << "----------------- slide ------------ " << std::endl;

      int sequence_length = sequence.length();
      int offset_max = sequence.length() - n_residues;

      for (int offset=0; offset<=offset_max; offset++) {
         int n_scored_residues = 0;
         double sum_score = 0;
         std::string running_sequence;
         for (int ires=0; ires<n_residues; ires++) {
            if ((ires+offset) < sequence_length) {
               mmdb::Residue *residue_p = scored_residues[ires].first;
               const std::map<std::string, double> &scored_map = scored_residues[ires].second;

               if (false) { // debug
                  std::map<std::string, double>::const_iterator it_debug;
                  for (it_debug=scored_map.begin();
                       it_debug!=scored_map.end();
                       it_debug++) {
                     std::cout << "   " << it_debug->first << " : " << it_debug->second << ", ";
                  }
                  std::cout << "\nDone map: " << std::endl;
               }

               char letter = sequence[ires+offset];
               std::string res_type = util::single_letter_to_3_letter_code(letter);
               // std::cout << "----------- debug:: res_type from " << letter << " is \"" << res_type << "\"" << std::endl;
               if (! res_type.empty()) {
                  std::map<std::string, double>::const_iterator it = scored_map.find(res_type);
                  if (it != scored_map.end()) {
                     const double &score = it->second;
                     running_sequence += letter;
                     sum_score += score;
                     n_scored_residues++;
                     if (false)
                        std::cout << "debug:: adding "
                                  << std::right << std::fixed << std::setw(9) << std::setprecision(4)
                                  << score << " for " << letter << " " << "ires " << ires << " "
                                  << residue_spec_t(residue_p) << std::endl;
                  } else {
                     if (false) { // debug
                        std::cout << "DEBUG:: Failed to find " << res_type << " in this map: " << std::endl;
                        std::map<std::string, double>::const_iterator it_debug;
                        std::cout << "This was the map: POINT B, it has size " << scored_map.size() << std::endl;
                        for (it_debug=scored_map.begin();
                             it_debug!=scored_map.end();
                             it_debug++) {
                           std::cout << "   " << it_debug->first << " : " << it_debug->second << ", ";
                        }
                        std::cout << "\nDone map: " << std::endl;
                     }
                     running_sequence += '.';
                  }
               } else {
                  running_sequence += '-';
               }
            }
         }
         if (n_scored_residues == n_residues) {
            results_t result(offset, sum_score, n_scored_residues, running_sequence, sequence_name,
                             true_sequence);
            results.push_back(result);
            std::cout << "INFO:: offset " << offset << " sum_score " << std::setw(8) << sum_score
                      << " n_scored_residues " << n_scored_residues << " " << running_sequence
                      << " gene-name " << gene_name
                      << " true-sequence " << true_sequence << std::endl;
         }
      }
   }
   get_results_addition_lock();
   results_container[sequence_name] = results;
   release_results_addition_lock();
}


coot::side_chain_densities::results_t
coot::side_chain_densities::get_result() const {

   double best_score = -9e10;
   results_t best_results;

   std::map<std::string, std::vector<results_t> >::const_iterator it;
   for (it=results_container.begin(); it!=results_container.end(); ++it) {
      const std::vector<results_t> &v = it->second;
      for (unsigned int i=0; i<v.size(); i++) {
         const results_t &result = v[i];
         if (result.sum_score > best_score) {
            best_score = result.sum_score;
            best_results = result;
         }
      }
   }

   return best_results;
}

bool
coot::side_chain_densities::test_grid_point_to_coords_interconversion() const {

   bool success = true;
   int n_per_side = 2 * n_steps + 1;
   for (int ix= -n_steps; ix<=n_steps; ix++) {
      for (int iy= -n_steps; iy<=n_steps; iy++) {
         for (int iz= -n_steps; iz<=n_steps; iz++) {
            int idx =
               (ix + n_steps) * n_per_side * n_per_side +
               (iy + n_steps) * n_per_side +
               (iz + n_steps);

            std::tuple<int, int, int> t = grid_index_to_grid(idx);
            if (std::get<0>(t) == ix)
               if (std::get<1>(t) == iy)
                  if (std::get<2>(t) == iz)
                     if (false) // noise - if the function works :-)
                        std::cout << "PASS " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                                  << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                                  << std::endl;

            if (std::get<0>(t) != ix) {
               std::cout << "FAIL " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                         << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                         << std::endl;
               success = false;
            }
            if (std::get<1>(t) != iy) {
               std::cout << "FAIL " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                         << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                         << std::endl;
               success = false;
            }
            if (std::get<2>(t) != iz) {
               std::cout << "FAIL " << std::setw(2) << ix << " " << std::setw(2) << iy << " " << std::setw(2) << iz << " idx " << idx
                         << " decode " << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t)
                         << std::endl;
               success = false;
            }
         }
      }
   }

   return success;

}


// gen_pts_file_name is optional argument for generating initial usable grid points.
// Calling function must clear the returned density_box_t
//
coot::density_box_t
coot::side_chain_densities::sample_map(mmdb::Residue *residue_this_p,
                                       mmdb::Residue *residue_next_p,
                                       mode_t mode,
                                       const clipper::Coord_orth &cb_pt,
                                       const std::vector<clipper::Coord_orth> &axes,
                                       const clipper::Xmap<float> &xmap,
                                       std::string gen_pts_file_name) const {


   bool gen_usable_points_flag = false;
   if (mode == GEN_USABLE_POINTS) gen_usable_points_flag = true;

   // 3 modes:
   //
   // (1) Generating the grid point list
   // (2) Sampling the maps to provide data for the distributions
   // (3) Sampling the user's map to generate a grid of data to test against the distributions.

   // This function is a bit ugly, because it used to generate grid points
   // and to sample the map - it is not split into 2 or 3 functions - which would
   // be cleaner, because I want only one version of the code that converts
   // from grid indices to 3D points.

   float step_size = grid_box_radius/static_cast<float>(n_steps);
   int n_per_side = 2 * n_steps + 1;
   int n_box_vol = n_per_side * n_per_side * n_per_side;

   if (! residue_this_p) return density_box_t(0,0,0);
   if (axes.empty())     return density_box_t(0,0,0);

   std::string res_name = residue_this_p->GetResName();
   std::string rot_name = get_rotamer_name(residue_this_p); // doesn't matter for SAMPLE_FOR_RESIDUE mode

   if (mode == SAMPLE_FOR_DB) {
      if (rot_name.empty()) {
         if (res_name == "GLY") {
            rot_name = "pseudo";
         } else {
            return density_box_t(0,0,0);
         }
      }
   }

   clipper::Coord_orth ca_pt(-1,-1,-1);

   int n_atoms = residue_this_p->GetNumberOfAtoms();
   std::vector<clipper::Coord_orth> residue_atom_positions;
   std::vector<std::pair<double, clipper::Coord_orth> > main_chain_atom_positions;
   residue_atom_positions.reserve(n_atoms);
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = residue_this_p->GetAtom(i);
      if (! at->isTer()) {
         clipper::Coord_orth pos = co(at);
         residue_atom_positions.push_back(pos);
         std::string atom_name(at->name);
         if (atom_name == " N  " || atom_name == " C  " || atom_name == " O  " ||
             atom_name == " H  " || atom_name == " CA ") {
            double r = 2.8;
            if (atom_name == " CA ")
               r = 1.6;
            std::pair<double, clipper::Coord_orth> p(r, pos);
            main_chain_atom_positions.push_back(p);
         }
         if (atom_name == " CA ") ca_pt = co(at);
      }
   }
   if (gen_usable_points_flag) {
      // need to reject points around N of next residue too
      n_atoms = residue_next_p->GetNumberOfAtoms();
      for (int i=0; i<n_atoms; i++) {
         mmdb::Atom *at = residue_next_p->GetAtom(i);
         if (! at->isTer()) {
            std::string atom_name(at->name);
            if (atom_name == " N  ") {
               clipper::Coord_orth pos = co(at);
               std::pair<double, clipper::Coord_orth> p(3.0, pos);
               main_chain_atom_positions.push_back(p);
            }
         }
      }
   }

   stats::single block_stats;

   std::ofstream f; // for generating the (index) side chain points file
   if (gen_usable_points_flag)
      f.open(gen_pts_file_name.c_str());

   // use shared pointer
   float *density_box = new float[n_box_vol];

   float unset_value = -1001.1;
   for (int i=0; i<n_box_vol; i++) density_box[i] = unset_value;

   // make main_chain_clashing_points - these should be the same for every residue
   // i.e. use something like this to make a reference list once
   // the grid size and n_steps is optimized.
   //
   // Note:  exclude grid points that are within 3(?)A of C, O, N and next N
   //
   std::set<int> main_chain_clashing_points;
   if (gen_usable_points_flag) {
      for (int ix= -n_steps; ix<=n_steps; ix++) {
         for (int iy= -n_steps; iy<=n_steps; iy++) {
            for (int iz= -n_steps; iz<=n_steps; iz++) {
               int idx =
                  (ix + n_steps) * n_per_side * n_per_side +
                  (iy + n_steps) * n_per_side +
                  (iz + n_steps);
               clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
               clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;
               clipper::Coord_orth ca_to_cb = cb_pt - ca_pt;
               clipper::Coord_orth ca_to_pt = pt_grid_point - ca_pt;
               double dp = clipper::Coord_orth::dot(ca_to_cb, ca_to_pt);
               if (dp < 0.0)
                  main_chain_clashing_points.insert(idx);
            }
         }
      }
   }

   for (int ix= -n_steps; ix<=n_steps; ix++) {
      for (int iy= -n_steps; iy<=n_steps; iy++) {
         for (int iz= -n_steps; iz<=n_steps; iz++) {
            int idx =
               (ix + n_steps) * n_per_side * n_per_side +
               (iy + n_steps) * n_per_side +
               (iz + n_steps);
            // note: when generating useable points, this set is empty
            if (gen_usable_points_flag || useable_grid_points.find(idx) != useable_grid_points.end()) {
               clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
               clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;
               if (gen_usable_points_flag) {
                  if (! is_close_to_atoms(main_chain_atom_positions, pt_grid_point)) {
                     if (in_sphere(pt_grid_point, cb_pt, grid_box_radius)) {
                        if (main_chain_clashing_points.find(idx) == main_chain_clashing_points.end()) {
                           if (gen_usable_points_flag) { // setting up for sampling
                              f << "setting grid point " << idx << " at "
                                << pt_grid_point.x() <<  " "
                                << pt_grid_point.y() <<  " "
                                << pt_grid_point.z() << std::endl;
                           }
                        }
                     }
                  }
               } else {
                  // more usual, the decision about acceptable grid points have already
                  // been made (some time ago).
                  float dv = util::density_at_point_by_linear_interpolation(xmap, pt_grid_point);
                  density_box[idx] = dv;
                  block_stats.add(dv);
                  if (false)
                     std::cout << "sample_map(): for block_stats " << pt_grid_point.format()
                               <<  " " << idx << " " << dv << "\n";
               }
            }
         }
      }
   }

   if (false)
      std::cout << "debug: block stats " << coot::residue_spec_t(residue_this_p)
                << " size " << block_stats.size()
                << " mean " << block_stats.mean() << " sd " << sqrt(block_stats.variance()) << std::endl;


   density_box_t db(density_box, residue_this_p, n_steps);

   // return negative values on failure
   std::tuple<double, double, double> ca_stats = get_stats_around_ca(residue_this_p, axes,
                                                                     0.5 * step_size, xmap);
   if (std::get<1>(ca_stats) > 0)
      db.set_around_ca_stats(std::get<0>(ca_stats), std::get<1>(ca_stats), std::get<2>(ca_stats));

   return db;
}

// return negative values on failure
std::tuple<double, double, double>
coot::side_chain_densities::get_stats_around_ca(mmdb::Residue *residue_this_p,
                                                const std::vector<clipper::Coord_orth> &axes,
                                                float step_size,
                                                const clipper::Xmap<float> &xmap) const {

   double mean = 0;
   double var = 0;
   double mean_of_positives = 0;

   int n_atoms = residue_this_p->GetNumberOfAtoms();
   mmdb::Atom *ca_at = 0;
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = residue_this_p->GetAtom(i);
      std::string atom_name(at->GetAtomName());
      if (atom_name == " CA ") {
         ca_at = at;
         break;
      }
   }
   if (ca_at) {

      stats::single s;
      stats::single s_above_zero;
      clipper::Coord_orth ca_pos = co(ca_at);
      int ilim = n_steps * n_steps;
      for (int ix= -n_steps; ix<=n_steps; ix++) {
         for (int iy= -n_steps; iy<=n_steps; iy++) {
            for (int iz= -n_steps; iz<=n_steps; iz++) {
               if ((ix * ix + iy * iy + iz * iz) <= ilim) {
                  clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
                  clipper::Coord_orth pt_grid_point = ca_pos + pt_in_grid;
                  float dv = util::density_at_point_by_linear_interpolation(xmap, pt_grid_point);
                  s.add(dv);
                  if (dv > 0.0)
                     s_above_zero.add(dv);
               }
            }
         }
      }

      mean = s.mean();
      var = s.variance();
      mean_of_positives = s_above_zero.mean();

   } else {
      // failure status
      mean = -1;
      var = -1;
      mean_of_positives = -1;
   }

   return std::tuple<double, double, double> (mean, var, mean_of_positives);
}


void
coot::side_chain_densities::normalize_density_boxes(const std::string &id) {

   // We are normalizing the boxes from our map/model (not reference data)

   // normalize_density_boxes_v2(id);

   normalize_density_boxes_v3(id); // Use CA density for scaling
}

void
coot::side_chain_densities::normalize_density_boxes_v1(const std::string &id) {

   // hacketty-hack scaling

   float sum = 0;
   float sum_sq = 0;
   int n_grid_pts = 0;

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      int nnn = density_boxes[i].nnn();
      const density_box_t &db = density_boxes[i];
      for (int j=0; j<nnn; j++) {
         if (db[j] > 0.0) {
            sum += db[j];
            sum_sq += db[j] * db[j];
            n_grid_pts++;
         }
      }
   }

   if (n_grid_pts > 0) {
      float mean = sum/static_cast<float>(n_grid_pts);
      float var = sum_sq/static_cast<float>(n_grid_pts) - mean * mean;
      float sd = sqrt(var);
      float scale_factor = 1.0/mean;
      std::cout << "Dataset from " << id << " mean " << mean << " scale_factor "
                << scale_factor << std::endl;
      for (std::size_t i=0; i<density_boxes.size(); i++) {
         density_box_t &db = density_boxes[i];
         db.scale_by(scale_factor); // don't scale "below zero" points
      }
   }
}

void
coot::side_chain_densities::normalize_density_boxes_v2(const std::string &id) {

   // Make the RMSd be 1.0 then the the mean to zero

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      int n_grid_pts = 0;
      float sum = 0;
      float sum_sq = 0;
      int nnn = density_boxes[i].nnn();
      density_box_t &db = density_boxes[i];
      for (int j=0; j<nnn; j++) {
         if (db[j] > -1000.0) {
            sum += db[j];
            sum_sq += db[j] * db[j];
            n_grid_pts++;
         }
      }

      if (n_grid_pts > 0) {
         float mean = sum/static_cast<float>(n_grid_pts);
         float var = sum_sq/static_cast<float>(n_grid_pts) - mean * mean;
         if (var > 0.0) {
            float sd = sqrt(var);
            float scale_factor = 1.0/sd;
            for (int j=0; j<nnn; j++) {
               if (db[j] > -1000.0) {
                  db.density_box[j] *= scale_factor;
               }
            }
            sum = 0;
            for (int j=0; j<nnn; j++) {
               if (db[j] > -1000.0) {
                  sum += db[j];
               }
            }
            mean = sum/static_cast<float>(n_grid_pts);
            for (int j=0; j<nnn; j++) {
               if (db[j] > -1000.0) {
                  db.density_box[j] -= mean;
               }
            }
         }
      }
   }
}

void
coot::side_chain_densities::normalize_density_boxes_v3(const std::string &id) {

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      int nnn = density_boxes[i].nnn();
      density_box_t &db = density_boxes[i];
      db.normalize_using_ca_stats();
   }
}

void
coot::density_box_t::normalize_using_ca_stats() {

   if (var_around_ca > 0) {
      int n = nnn();
      if (mean_of_positives_around_ca > 0.0) {
         if (false)
            std::cout << "Debug: normalizing with mean of CA positives "
                      << mean_of_positives_around_ca << " for "
                      << residue_spec_t(residue_p) << std::endl;
         float sf = 0.995 * 1.0/mean_of_positives_around_ca;
         unsigned int n_scaled = 0;
         for (int j=0; j<n; j++) {
            if (density_box[j] > -1000.0) {
               density_box[j] *= sf;
               n_scaled++;
            }
         }
         // std::cout << "DEBUG:: n_scaled " << n_scaled << std::endl;
      } else {
         // mark as a baddie - perhaps this
         // should have its own flag?
         var_around_ca = -1.0;
         is_weird = true;
      }
   } else {
      std::string rn;
      is_weird = true;
      if (residue_p)
         rn = residue_p->GetResName();
      std::cout << "ERROR:: Failed variance in normalize_using_ca_stats() "
                << residue_spec_t(residue_p) << " " << rn << std::endl;
   }

}



void
coot::side_chain_densities::write_density_boxes() const {

   for (std::size_t i=0; i<density_boxes.size(); i++) {
      const density_box_t &db = density_boxes[i];
      // don't write boxes with weird/non-filled stats
      if (db.var_around_ca > 0.0)
         if (! db.is_weird)
            write_density_box(db, id);
   }
}

std::pair<float, float>
coot::density_box_t::mean_and_variance() const {

   // of non-masked points above the zero

   stats::single s;
   int n = 2 * n_steps + 1;
   int count = 0;
   int nnn = n * n * n;
   for (int i=0; i<nnn; i++) {
      const float &d = density_box[i];
      if (d > 0.0) {
         s.add(d);
         count++;
      }
   }

   float mean = 99999999999.9;
   float var = 0;
   if (count > 0) {
      mean = s.mean();
      var = s.variance();
   }
   return std::pair<float, float> (mean, var);
}

void
coot::side_chain_densities::write_density_box(const coot::density_box_t &db, const std::string &id) const {

   float *density_box = db.density_box;
   int n_steps = db.n_steps;
   mmdb::Residue *residue_p = db.residue_p;

   if (! residue_p) return;
   std::string res_name = residue_p->GetResName();
   std::string rot_name = get_rotamer_name(residue_p);
   std::string dir = "side-chain-data";

   if (!rot_name.empty()) {
      std::string rot_dir = dir + "/" + res_name + "/" + rot_name;
      std::string file_name = rot_dir + "/" + id + "-" + residue_p->GetChainID() + util::int_to_string(residue_p->GetSeqNum()) + ".tab";

      if (! file_exists(rot_dir))
         util::create_directory(rot_dir);
      // std::cout << "write_density_box() to filename " << file_name << std::endl;

      std::ofstream f(file_name.c_str());
      if (f) {
         int n_per_side = 2 * n_steps + 1;
         int n_box_vol = n_per_side * n_per_side * n_per_side;
         for (int i=0; i<n_box_vol; i++) {
            float v = density_box[i];
            if (clipper::Util::is_nan(v)) {
               std::cout << "ERROR:: " << file_name << " found a nan " << i << std::endl;
               exit(1);
            } else {
               f << v << " ";
               if (i%n_per_side == 0)
                  f << "\n";
            }
         }
         f << "\n";
      } else {
         std::cout << "WARNING:: cannot open file " << file_name << std::endl;
      }
   } else {
      // these need investigation - Coot bug?
      bool verbose = false;
      if (verbose)
         std::cout << "WARNING:: no rotamer name for " << id << " " << residue_spec_t(residue_p)
                   << " " << res_name << std::endl;
   }
}

// stores density boxes - change the name of the function
void
coot::side_chain_densities::proc_chain(const std::string &id, mmdb::Chain *chain_p,
                                       const clipper::Xmap<float> &xmap) {

   int n_residues = chain_p->GetNumberOfResidues();
   int last_res = n_residues - 2;
   for (int ires=1; ires<=last_res; ires++) {
      mmdb::Residue *t = chain_p->GetResidue(ires);
      if (t) {

        std::string rn = t->GetResName();
        if (rn == "UNK") continue;

        if (! util::is_standard_amino_acid_name(rn)) continue;

         // don't forget that ALA are useful to search, as is GLY, but
         // that will need a special function to find an imaginary CB position
         //
         std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
            get_residue_axes(t);
         const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
         if (! axes.empty()) {
            // sample_masked density around 8A around CB
            clipper::Coord_orth cb_pt = cb_pos_and_axes.first;
            bool gen_flag = false;
            mode_t mode = SAMPLE_FOR_DB;
            density_box_t db = sample_map(t, 0, mode, cb_pt, axes, xmap);
            if (! db.empty()) {
               // std::cout << "Storing density box for residue " << residue_spec_t(t) << " "
               //           << residue_spec_t(db.residue_p) << std::endl;
               store_density_box(db); // push back to density_boxes vector
            }
         }
      }
   }
}

#include "coot-utils/c-beta-deviations.hh"

std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> >
coot::side_chain_densities::get_residue_axes_type_GLY(mmdb::Residue *this_residue) const {

   mmdb::Atom *CA_at = 0;
   mmdb::Atom *N_at  = 0;
   mmdb::Atom *C_at  = 0;
   int n_atoms = this_residue->GetNumberOfAtoms();
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = this_residue->GetAtom(i);
      std::string atom_name = at->name;
      std::string alt_loc = at->altLoc;
      if (! at->isTer()) {
         if (alt_loc.empty()) {
            if (atom_name == " CA ") CA_at = at;
            if (atom_name == " N  ")  N_at = at;
            if (atom_name == " C  ")  C_at = at;
         }
      }
   }
   if (N_at && CA_at && C_at) {
      atom_quad q(N_at, CA_at, C_at, 0); // only 3 atoms used
      clipper::Coord_orth cb_pos = make_CB_ideal_pos(q, "ALA");
      clipper::Coord_orth ca_pos = co(CA_at);
      clipper::Coord_orth c_pos = co(C_at);
      clipper::Coord_orth n_pos = co(N_at);
      std::vector<clipper::Coord_orth> axes = make_axes(ca_pos, cb_pos, c_pos, n_pos);
      return std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > (cb_pos, axes);
   } else {
      std::cout << "ERROR:: BAD GLY " << residue_spec_t(this_residue) << std::endl;
      std::vector<clipper::Coord_orth> axes;
      clipper::Coord_orth cb_pos(0,0,0);
      return std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > (cb_pos, axes);
   }
   
}

std::vector<clipper::Coord_orth>
coot::side_chain_densities::make_axes(const clipper::Coord_orth &pt_ca_this,
                                      const clipper::Coord_orth &pt_cb_this,
                                      const clipper::Coord_orth &pt_c_this,
                                      const clipper::Coord_orth &pt_n_this) const {

   std::vector<clipper::Coord_orth> v;
   // Here add a bit of noise to pt_cb_this to "sample"
   // density from non-ideal orientation

   clipper::Coord_orth axis_1((pt_cb_this - pt_ca_this).unit());
   clipper::Coord_orth nc(pt_c_this - pt_n_this);
   clipper::Coord_orth nc_uv(nc.unit());
   clipper::Coord_orth cp_1(clipper::Coord_orth::cross(axis_1, nc_uv));
   clipper::Coord_orth cp_2(clipper::Coord_orth::cross(cp_1, axis_1));

   clipper::Coord_orth axis_2(cp_1.unit());
   clipper::Coord_orth axis_3(cp_2.unit());

   double dp = clipper::Coord_orth::dot(nc_uv, axis_3);
   double theta = acos(dp);
   // std::cout << residue_spec_t(this_residue) << " "
   // << clipper::Util::rad2d(theta) << std::endl;
   v.push_back(axis_1);
   v.push_back(axis_2);
   v.push_back(axis_3);

   return v;
}


std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> >
coot::side_chain_densities::get_residue_axes(mmdb::Residue *residue_p) const {

   // only look for atoms with alt conf ""

   std::string res_name(residue_p->GetResName());
   if (res_name == "GLY") return get_residue_axes_type_GLY(residue_p);

   // OK, we have a CB (presumably)

   std::vector<clipper::Coord_orth> v;
   clipper::Coord_orth cb_pos(0,0,0);

   mmdb::Atom *CA_at = 0;
   mmdb::Atom *CB_at = 0;
   mmdb::Atom *N_at  = 0;
   mmdb::Atom *C_at  = 0;
   mmdb::Atom *O_at  = 0;

   int n_atoms = residue_p->GetNumberOfAtoms();
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = residue_p->GetAtom(i);
      std::string atom_name = at->name;
      std::string alt_loc = at->altLoc;
      if (! at->isTer()) {
         if (alt_loc.empty()) {
            if (atom_name == " CA ") CA_at = at;
            if (atom_name == " CB ") CB_at = at;
            if (atom_name == " N  ")  N_at = at;
            if (atom_name == " C  ")  C_at = at;
            if (atom_name == " O  ")  O_at = at;
         }
      }
   }

   if (CA_at && CB_at && N_at && C_at && O_at) {

         clipper::Coord_orth pt_ca_this = co(CA_at);
         clipper::Coord_orth pt_cb_this = co(CB_at);
         clipper::Coord_orth pt_c_this = co(C_at);
         clipper::Coord_orth pt_n_this = co(N_at);

         v = make_axes(pt_ca_this, pt_cb_this, pt_c_this, pt_n_this);
         cb_pos = co(CB_at);

   }
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > p(cb_pos, v);
   return p;
}

void
coot::side_chain_densities::gen_useable_grid_points(mmdb::Residue *residue_this_p,
                                                    mmdb::Residue *residue_next_p,
                                                    int n_steps, float grid_box_radius,
                                                    const std::string &gen_pts_file_name) const {

   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_this_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   if (! axes.empty()) {
      // sample_masked density around 8A around CB
      clipper::Coord_orth cb_pt = cb_pos_and_axes.first;

      // make a block
      mode_t mode = GEN_USABLE_POINTS;
      clipper::Xmap<float> dummy;
      density_box_t block = sample_map(residue_this_p, residue_next_p, mode, cb_pt,
                                       axes, dummy, gen_pts_file_name);
   }
}

std::pair<std::string, std::string>
coot::side_chain_densities::map_key_to_residue_and_rotamer_names(const std::string &key) const {

   std::string::size_type pos = key.find_last_of(":");
   std::string residue_name = key.substr(0, pos);
   std::string rotamer_name = key.substr(pos+1, key.length());

   return std::pair<std::string, std::string>(residue_name, rotamer_name);
}

void
coot::side_chain_densities::fill_residue_blocks(mmdb::Manager *mol, const std::string &chain_id,
                                                int resno_start, int resno_end,
                                                const clipper::Xmap<float> &xmap) {

   const std::vector<mmdb::Residue *> residues = make_a_run_of_residues(mol, chain_id, resno_start, resno_end);
   fill_residue_blocks(residues, xmap);

}

// we can do better normalization of the grids for the user/test structure if
// we do them all at the same time.
void
coot::side_chain_densities::fill_residue_blocks(const std::vector<mmdb::Residue *> &residues,
                                                const clipper::Xmap<float> &xmap) {

   if (density_block_map_cache.size() > 0) {
      return; // already done
   } else {
      get_results_addition_lock();
      for (std::size_t i=0; i<residues.size(); i++) {
         mmdb::Residue *residue_p = residues[i];
         mode_t mode = SAMPLE_FOR_RESIDUE;
         std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
            get_residue_axes(residue_p);
         const clipper::Coord_orth &cb_pt = cb_pos_and_axes.first;
         const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
         density_box_t block = sample_map(residue_p, 0, mode, cb_pt, axes, xmap);
         block.normalize_using_ca_stats();
         density_block_map_cache[residue_p] = block;
      }
      add_mean_and_variance_to_individual_density_blocks();
      release_results_addition_lock();
   }

}

void
coot::side_chain_densities::normalize_density_blocks() {

   std::map<mmdb::Residue *, density_box_t>::const_iterator it;
   unsigned int n_grid_pts = 0;
   double sum = 0;

   stats::single s_all;
   for(it=density_block_map_cache.begin(); it!=density_block_map_cache.end(); it++) {
      const density_box_t &block = it->second;
      if (! block.empty()) {
         int nnn = block.nnn();
         for (int i=0; i<nnn; i++) {
            if (block[i] > 0.0) {
               sum += block[i];
               n_grid_pts++;
            }
            s_all.add(block[i]);
         }
      }
   }
   if (n_grid_pts > 0) {
      double av = sum/static_cast<double>(n_grid_pts);
      double sc = mn_scale_for_normalized_density/av;
      std::map<mmdb::Residue *, density_box_t>::iterator it_inner;
      for(it_inner=density_block_map_cache.begin(); it_inner!=density_block_map_cache.end(); it_inner++) {
         density_box_t &block = it_inner->second; // because not const
         block.scale_by(sc);
      }
   }
}

void
coot::side_chain_densities::add_mean_and_variance_to_individual_density_blocks() {

   stats::single s;
   stats::single s_positive;
   std::map<mmdb::Residue *, density_box_t>::iterator it;
   for(it=density_block_map_cache.begin(); it!=density_block_map_cache.end(); it++) {
      density_box_t &block = it->second;
      if (! block.empty()) {
         int nnn = block.nnn();
         for (int i=0; i<nnn; i++) {
            const float &bi = block[i];
            if (bi > -1000.0)
               s.add(bi);
            if (bi > 0.0)
               s_positive.add(bi);
         }
         if (false)
            std::cout << "debug:: in add_mean_and_variance_to_individual_density_blocks() "
                      << residue_spec_t(it->first) << " mean: " << s.mean() << " variance: " << s.variance()
                      << std::endl;
         block.set_stats(s.mean(), s.variance(), s_positive.mean());
      }
   }
}

coot::density_box_t
coot::side_chain_densities::get_block(mmdb::Residue *residue_p) const {

   std::map<mmdb::Residue *, density_box_t>::const_iterator it;
   it = density_block_map_cache.find(residue_p); // this cannot (must not) fail - so make
                                                 // sure that fill_residue_blocks is called
                                                 // before this function is called.
   if (it == density_block_map_cache.end()) {
      std::cout << "ERROR:: in get_block(): Hideous failure!" << std::endl;
   }
   return it->second;
}



// the given residue needs to have a CB - caller should check and make one if
// the model doesn't have one
// limit_to_correct_rotamers_only is optional arg, default false, and is for debugging the llr for
// correct solutions - which is bad/low and why?
// verbose_output_mode is optional arg, default true
//
std::map<std::string, double>
coot::side_chain_densities::get_rotamer_likelihoods(mmdb::Residue *residue_p,
                                                    const clipper::Xmap<float> &xmap,
                                                    bool limit_to_correct_rotamers_only,
                                                    bool verbose_output_mode) {

   // To see what's going on with residue scoring.
   // verbose_output_mode = true;

   // fill_residue_blocks() has been called before we get here

   std::map<std::string, double> bs; // return this, best_score_for_res_type

   if (density_block_map_cache.size() == 0) {
      std::cout << "ERROR:: Cache is empty - fill it first" << std::endl;
      return bs;
   }

   // are axes needed?
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   if (! axes.empty()) {
      // sample_masked density around CB

      // retrieve the block for this residue (filled in fill_residue_blocks())

      density_box_t block = get_block(residue_p);

      if (verbose_output_mode)
         std::cout << "debug:: in get_rotamer_likelihoods() residue: "
                   << residue_spec_t(residue_p) << " block has mean " << block.mean
                   << " and var " << block.var << " sd " << sqrt(block.var) << std::endl;

      if (block.empty()) {

         std::cout << "WARNING:: failed to get a density block for "
                   << residue_spec_t(residue_p) << std::endl;

      } else {

         if (block.is_weird) {
            // do nothing
            std::cout << "WARNING:: weird density block for "
                      << residue_spec_t(residue_p) << std::endl;
         } else {

            // Happy Path

            // compare the block of density to every rotamer (of every residue)

            std::string res_name = residue_p->GetResName();
            std::string rot_name = get_rotamer_name(residue_p);

            if (verbose_output_mode)
               std::cout << "debug:: in get_rotamer_likelihoods(): residue "
                         << residue_spec_t(residue_p) << " " << res_name << " " << rot_name << std::endl;

            std::pair<bool, std::vector<std::pair<std::string, std::string> > > rotamer_limits;
            rotamer_limits.first = false;

            // just test the rotamer that this residue is, combine with
            // outputting the likelihoods for the grids in get_log_likelihood_ratio()
            // in the "check the engine" ouput.

            if (limit_to_correct_rotamers_only)
               rotamer_limits.first = true; // testing the likelihoods of the correct rotamers

            rotamer_limits.second.push_back(std::pair<std::string, std::string> (res_name, rot_name));

            std::map<std::string, double> probability_map =
               compare_block_vs_all_rotamers(block, residue_p, data_dir, rotamer_limits, xmap);

            // find the max value so that I can mark it and close others for screen output
            double best_score = -999999999999999.9;
            std::map<std::string, double>::const_iterator it;
            for (it=probability_map.begin(); it!=probability_map.end(); it++) {
               if (it->second > best_score) {
                  best_score = it->second;
               }
            }

            std::map<std::string, double> best_score_for_res_type;
            // std::cout << "here with probability_map size " << probability_map.size() << std::endl;
            for (it=probability_map.begin(); it!=probability_map.end(); it++) {
               const std::string &key = it->first; // resname:rot_name
               // std::cout << "key: " << key << std::endl;
               std::string m;
               std::pair<std::string, std::string> rrp = map_key_to_residue_and_rotamer_names(key);
               const std::string &res_name = rrp.first;
               const std::string &rot_name = rrp.second;
               const double &score = it->second;
               if (best_score_for_res_type.find(res_name) == best_score_for_res_type.end()) {
                  best_score_for_res_type[res_name] = score;
                  if (false)
                     std::cout << "DEBUG:: first  score for res-type " << std::setw(10) << std::left << key << " "
                               << std::fixed << std::right << std::setw(7) << std::setprecision(2)
                               << score << std::endl;
               } else {
                  if (score > best_score_for_res_type[res_name]) {
                     best_score_for_res_type[res_name] = score;
                     if (false)
                        std::cout << "DEBUG:: better score for res-type " << std::setw(10) << std::left << key << " "
                                  << std::fixed << std::right << std::setw(7) << std::setprecision(2)
                                  << score << std::endl;
                  }
               }
            }
            bs = best_score_for_res_type;

            if (verbose_output_mode) {
               for (it=best_score_for_res_type.begin(); it!=best_score_for_res_type.end(); it++) {
                  const std::string &res_type = it->first;
                  const double &score = it->second;
                  std::string m;
                  if (best_score < 0) {
                     if (score > 1.1 * best_score) m = " ooo";
                  } else {
                     if (score > 0)
                        if (score > 0.5 * best_score)
                           m = " ooo";
                  }
                  if (score == best_score) m = " ***";
                  std::cout << "   " << res_type << " "
                            << std::fixed << std::right << std::setw(8) << std::setprecision(2)
                            << score << m << std::endl;
               }
            }
         }
      }
   }

   return bs;
}


std::map<std::string, double>
coot::side_chain_densities::compare_block_vs_all_rotamers(density_box_t block,
                                                          mmdb::Residue *residue_p, // for debugging
                                                          const std::string &data_dir,
                                                          const std::pair<bool, std::vector<std::pair<std::string, std::string> > > &rotamer_limits,
                                                          const clipper::Xmap<float> &xmap) {

   std::map<std::string, double> probability_map;

   std::string glob_pattern = "*";
   std::vector<std::string> dirs = coot::util::glob_files(data_dir, glob_pattern);
   // std::cout << "found " << dirs.size() << " files in " << data_dir << std::endl;

   for (std::size_t idir=0; idir<dirs.size(); idir++) {
      const std::string &res_dir = dirs[idir];

      std::vector<std::string> rot_dirs = coot::util::glob_files(res_dir, glob_pattern);
      for (std::size_t jdir=0; jdir<rot_dirs.size(); jdir++) {
         const std::string &rot_dir = rot_dirs[jdir];

         std::string res = util::file_name_non_directory(res_dir);
         std::string rot = util::file_name_non_directory(rot_dir);
         std::string key = res + ":" + rot;

         bool do_it = false;

         if (false)
            std::cout << "debug:: in compare_block_vs_all_rotamers() rotamer limits first "
                      << rotamer_limits.first << std::endl;
         if (! rotamer_limits.first) {
            // so, don't apply the limits
            do_it = true;
         } else {
            if (rotamer_limits.first) {
               // testing path
               if (! rotamer_limits.second.empty()) {
                  bool found = false;
                  for (std::size_t i=0; i<rotamer_limits.second.size(); i++) {
                     std::string limit_key = rotamer_limits.second[i].first + ":" + rotamer_limits.second[i].second;
                     if (limit_key == key) {
                        found = true;
                        break;
                     }
                  }
                  if (found)
                     do_it = true;
               }
            }
         }

         if (do_it) {

            if (false)
               std::cout << "debug:: in compare_block_vs_all_rotamers(): rot_dir: " << rot_dir << std::endl;
            // std::cout << "debug:: in compare_block_vs_all_rotamers(): block var " << block.var << std::endl;
            std::pair<bool, double> p = compare_block_vs_rotamer(block, residue_p, rot_dir, xmap);
            if (p.first) {
               // std::cout << "debug:: in compare_block_vs_rotamer() pr: " << key << " " << p.second
               //           << std::endl;

               probability_map[key] = p.second;
            }
         } else {
            std::cout << "debug:: in compare_block_vs_rotamer() do_it was false " << rot_dir << std::endl;
         }
      }
   }
   return probability_map;
}

bool
coot::side_chain_densities::get_test_map_is_above_model_mean(const unsigned int &grid_idx,
                                                             const density_box_t &block,
                                                             const double &mean) const {
   double x = block[grid_idx];
   return (x > mean);

}

double
coot::side_chain_densities::get_log_likelihood(const unsigned int &grid_idx,
                                               const density_box_t &block,
                                               const double &mean,
                                               const double &variance,
                                               const double &skew) const {

   double x = block[grid_idx];
   double a = x - mean;
   double c_part = log(sqrt(1.0/(2.0 * M_PI * variance)));
   double e_part = -0.5*a*a/variance;

   if (false)
      std::cout << "get_log_likelihood() x " << x
                << " grid-idx " << grid_idx << " nz " << a/sqrt(variance) << std::endl;

   return c_part + e_part;
}

double
coot::side_chain_densities::get_grid_point_distance_from_grid_centre(const unsigned int &idx,
                                                                     const double &step_size) const {

   std::tuple<int, int, int> grid_coord = grid_index_to_grid(idx);
   // this is not the same axes system as the grid points - but they are on the same scale
   clipper::Coord_orth pt_in_grid(std::get<0>(grid_coord) * step_size,
                                  std::get<1>(grid_coord) * step_size,
                                  std::get<2>(grid_coord) * step_size);
   double l = std::sqrt(pt_in_grid.lengthsq());
   return l;
}

// log_likelihood ratio vs the gaussian sphere null hypothesis
double
coot::side_chain_densities::get_log_likelihood_ratio(const unsigned int &grid_idx,
                                                     const density_box_t &block,
                                                     const double &step_size,
                                                     const double &mean,
                                                     const double &variance_in,
                                                     const double &skew) const {

   double density_val = block[grid_idx];
   if (density_val > mn_density_block_sample_x_max)
      density_val = mn_density_block_sample_x_max;

   double variance = variance_in; // test/hack
   variance = 0.11; // observed variances are not useful.

   double var_scale = variance/block.var;
   double sd_scale = sqrt(var_scale);
   double mean_offset = mean - block.mean;

   if (false) {
      std::cout << "debug:: variance " << variance << " block.var " << block.var
                << " sd_scale " << sd_scale << std::endl;
      std::cout << "debug:: scaling density_val " << density_val << " with " << sd_scale
                << " and mean offset " << mean_offset << std::endl;
   }
   // double x = density_val * sd_scale - mean_offset;
   double x = density_val * 1.0; // weird scaling makes LLR better.

   double z = x - mean;
   double c_part = log(sqrt(1.0/(2.0 * M_PI * variance)));
   double e_part = -0.5*z*z/variance;

   // null hypothesis

   double nhs = null_hypothesis_scale;
   nhs = 2.0; // why is this 2.0?

   // distance between grid point and the CB
   double d = get_grid_point_distance_from_grid_centre(grid_idx, step_size);
   double x0_null_hypothesis = d;
   double c_part_null_normal = 1.0/(sqrt(2.0 * M_PI * null_hypothesis_sigma * null_hypothesis_sigma));
   double z_null_normal = d;
   double e_part_null_normal = exp(-((z_null_normal*z_null_normal)/(2.0 * null_hypothesis_sigma * null_hypothesis_sigma)));
   double x0_fake_density = nhs * c_part_null_normal * e_part_null_normal;

   double z_null = x0_fake_density - mean; // z value for the null hypothesis density value
   double e_part_normal = -0.5 * z_null * z_null / variance;

   if (false) {
      std::cout << "get_log_likelihood_ratio() x " << x
                << " grid-idx " << grid_idx << " nz " << z/sqrt(variance) << std::endl;
      std::cout << "in get_log_likelihood_ratio() null-hyp scale sigma " << null_hypothesis_scale
                << " " << null_hypothesis_sigma << std::endl;
      std::cout << "in get_log_likelihood_ratio() d " << d << std::endl;
      std::cout << "in get_log_likelihood_ratio() A " << e_part << " " << e_part_normal << std::endl;
      std::cout << "in get_log_likelihood_ratio() B " << z_null << " " << x0_fake_density << std::endl;
      std::cout << "in get_log_likelihood_ratio() C " << c_part_null_normal << " " << e_part_null_normal
                << std::endl;
   }
   double w = 1.0;
   // w = 2.3/(d + 1.0);
   // w = 1.0 - d * 0.166;
   // w = w * w * w;
   // w = 1.0;
   double diff = w * (e_part - e_part_normal);

   // remove this hideous baddies: Magic number - needs optimizing
   double mn_log_likelihood_ratio_difference_max = 18.0;
   if (diff < mn_log_likelihood_ratio_difference_min)
      diff = mn_log_likelihood_ratio_difference_min;
   if (diff > mn_log_likelihood_ratio_difference_max)
      diff = mn_log_likelihood_ratio_difference_max;

   if (false) // debug/check the engine
       /*
                << " e_part: " << std::setw(10) << e_part
                << " e_part_normal: " << std::setw(8) << e_part_normal
      */
      std::cout << "engine: idx: " << grid_idx
                << " with density_val " << std::setw(8) << std::right << std::setprecision(5) << density_val
                << " dv-mean " << std::right << std::setw(8) << std::setprecision(5) << z
                << " x0_fake_density " << std::setw(8) << x0_fake_density
                << " mean " << std::setw(5) << mean << " sigma " << sqrt(variance)
                << " return " << std::fixed << std::right << std::setprecision(6) << diff
                << "\n";

   return diff;
}

std::tuple<int, int, int>
coot::side_chain_densities::grid_index_to_grid(int grid_idx) const {

   int n_z = 0;
   int n_y = 0;
   int n_x = 0;
   int n = 2 * n_steps + 1;

   while (grid_idx >= (n * n)) {
      grid_idx -= n * n;
      n_x++;
   }
   while (grid_idx >= n) {
      grid_idx -= n;
      n_y++;
   }
   n_z = grid_idx;

   std::tuple<int, int, int> t(n_x - n_steps , n_y - n_steps, n_z - n_steps);
   return t;

}

// Manhattan test - not very useful
bool
coot::side_chain_densities::in_sphere(int grid_idx, const int &n_steps) const {

   bool inside = true;

   std::tuple<int, int, int> t = grid_index_to_grid(grid_idx);
   int n_x = std::get<0>(t);
   int n_y = std::get<1>(t);
   int n_z = std::get<2>(t);
   int delta = abs(n_x) + abs(n_y) + abs(n_z);
   if (delta > n_steps)
      inside = false;

   return inside;

}

bool
coot::side_chain_densities::in_sphere(const clipper::Coord_orth &pt,
                                      const clipper::Coord_orth &cb,
                                      const double &d_max) const {
   return ((pt-cb).lengthsq() < (d_max * d_max));
}

// We fill the rotamer grid cache, so not const.
// Maybe there should be 2 separate functions
//
std::pair<bool, double>
coot::side_chain_densities::compare_block_vs_rotamer(density_box_t block,
                                                     mmdb::Residue *residue_p,
                                                     const std::string &rotamer_dir,
                                                     const clipper::Xmap<float> &xmap) {

   // This function (-is-) was slow - conversion of strings to numbers
   // So, to fix that, in the calling function (compare_block_vs_all_rotamers()),
   // read in all the stats files - write and read the stats file as binaries
   // if it's still slow, so that this function is not needed.

   // I want to match the density mean and sigma of this blocks mean and sigma with the
   // overall stats with a scale (to match sigmas) and offset (to match means)

   if (false)
      std::cout << "debug:: compare_block_vs_rotamer() residue "
                << coot::residue_spec_t(residue_p) << " " << rotamer_dir
                << std::endl;

   bool do_debug_scoring = false; // count the number of times the the sample map is more than
                                  // the template stats mean. Should be about 50-50?
                                  // Turns out not. 100-150 seems to work.
   
   unsigned int n_grid_pts = 0;
   unsigned int n_above = 0;

   bool success = false; // initially
   double sum_log_likelihood = 0.0;
   double step_size = grid_box_radius/static_cast<float>(n_steps);

   std::map<std::string, std::map<unsigned int, std::tuple<double, double, double> > >::const_iterator it = rotamer_dir_grid_stats_map_cache.find(rotamer_dir);

   // std::cout << "------- calling get_log_likelihood_ratio() for rotamer_dir " << rotamer_dir << std::endl;

   if (it != rotamer_dir_grid_stats_map_cache.end()) {
      success = true;
      const std::map<unsigned int, std::tuple<double, double, double> > &stats_map = it->second;
      std::map<unsigned int, std::tuple<double, double, double> >::const_iterator it_inner;
      for (it_inner=stats_map.begin(); it_inner!=stats_map.end(); it_inner++) {
         const unsigned int &grid_idx = it_inner->first;
         const std::tuple<double, double, double> &m_v_s = it_inner->second;
         const double &mean = std::get<0>(m_v_s);
         const double &var  = std::get<1>(m_v_s);
         const double &skew = std::get<2>(m_v_s);
         // std::cout << "debug:: compare_block_vs_rotamer() A var " << var << " block.var " << block.var << std::endl;
         // double ll = get_log_likelihood(grid_idx, block, mean, var, skew);
         double llr = get_log_likelihood_ratio(grid_idx, block, step_size, mean, var, skew);
         sum_log_likelihood += llr;
         if (do_debug_scoring) {
            bool above = get_test_map_is_above_model_mean(grid_idx, block, mean);
            if (above)
               n_above++;
            n_grid_pts++;
         }
      }
   } else {
      std::string glob_pattern = "stats.table";
      std::vector<std::string> tables = coot::util::glob_files(rotamer_dir, glob_pattern);
      if (tables.size() == 1) {
         std::map<unsigned int, std::tuple<double, double, double> > stats_map;
         std::string stats_table_file_name = tables[0];

         // std::cout << "stats_table_file_name: " << stats_table_file_name << std::endl;
         std::ifstream f(stats_table_file_name.c_str());
         if (f) {
            std::string line;
            unsigned int n_grid_points = 0;
            while (std::getline(f, line)) {
               std::vector<std::string> words = coot::util::split_string_no_blanks(line);
               if (words.size() == 4) { // 5 with kurtosis
                  unsigned int grid_idx = util::string_to_int(words[0]);
                  double mean = util::string_to_double(words[1]);
                  double var  = util::string_to_double(words[2]);
                  double skew = util::string_to_double(words[3]);
                  if (false)
                     std::cout << "debug:: compare_block_vs_rotamer() B var " << var << " block.var " << block.var
                               << std::endl;
                  if (var < 0.0) std::cout << "ERROR:: negative variance " << var << std::endl;
                  // double ll = get_log_likelihood(grid_idx, block, mean, var, skew);
                  double llr = get_log_likelihood_ratio(grid_idx, block, step_size, mean, var, skew);
                  sum_log_likelihood += llr;
                  n_grid_points++;
                  if (do_debug_scoring) {
                     bool above = get_test_map_is_above_model_mean(grid_idx, block, mean);
                     if (above)
                        n_above++;
                     n_grid_pts++;
                  }
                  if (false)
                     std::cout << "for rotamer_dir " << rotamer_dir << " grid point "
                               << grid_idx << " density " << block[grid_idx]
                               << " adding " << llr << std::endl;

                  // cache that
                  std::tuple<double, double, double> t(mean, var, skew);
                  stats_map[grid_idx] = t;

               }
            }

            get_results_addition_lock(); // not really results...
            rotamer_dir_grid_stats_map_cache[rotamer_dir] = stats_map;
            release_results_addition_lock();

            if (n_grid_points > 0) {

               success = true;

               if (false)
                  std::cout << "debug:: " << std::left << std::setw(27) << rotamer_dir
                            << " sum-ll: " << std::right << sum_log_likelihood
                            << std::endl;
            } else {
               std::cout << "zero grid points " << rotamer_dir << std::endl;
            }

         } else {
            std::cout << "ERROR:: failed to open stats file " << stats_table_file_name
                      << std::endl;
         }
      }
   }

   if (do_debug_scoring) {
      std::cout << "debug:: above the mean " << n_above << " of " << n_grid_pts << std::endl;
   }

   return std::pair<bool, double>(success, sum_log_likelihood);
}

void
coot::side_chain_densities::check_useable_grid_points(mmdb::Residue *residue_p,
                                                      const std::string &useable_grid_points_mapped_to_residue_file_name) const {

   int n_per_side = n_steps * 2 + 1;
   float step_size = grid_box_radius/static_cast<float>(n_steps);
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   const clipper::Coord_orth &cb_pt = cb_pos_and_axes.first;

   std::ofstream f(useable_grid_points_mapped_to_residue_file_name.c_str());
   if (f) {
      if (! axes.empty()) {
         for (int ix= -n_steps; ix<=n_steps; ix++) {
            for (int iy= -n_steps; iy<=n_steps; iy++) {
               for (int iz= -n_steps; iz<=n_steps; iz++) {
                  int idx =
                     (ix + n_steps) * n_per_side * n_per_side +
                     (iy + n_steps) * n_per_side +
                     (iz + n_steps);
                  if (useable_grid_points.find(idx) != useable_grid_points.end()) {
                     clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
                     clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;
                     if (useable_grid_points.find(idx) != useable_grid_points.end()) {
                        f << "check-useable-grid-points " << idx << " "
                          << pt_grid_point.x() << " "
                          << pt_grid_point.y() << " "
                          << pt_grid_point.z() << "\n";
                     }
                  }
               }
            }
         }
      }
   }
   f.close();
}


void
coot::side_chain_densities::check_stats(mmdb::Residue *residue_p,
                                        const std::string &res_name,
                                        const std::string &rot_name) const {

   int n_per_side = n_steps * 2 + 1;
   float step_size = grid_box_radius/static_cast<float>(n_steps);
   std::pair<clipper::Coord_orth, std::vector<clipper::Coord_orth> > cb_pos_and_axes =
      get_residue_axes(residue_p);
   const std::vector<clipper::Coord_orth> &axes = cb_pos_and_axes.second;
   const clipper::Coord_orth &cb_pt = cb_pos_and_axes.first;
   if (! axes.empty()) {
      std::string dir = "side-chain-data"; // here use a class data item (which is set on construction)
      std::string rot_dir = dir + "/" + res_name + "/" + rot_name;
      std::string file_name = rot_dir + "/" + "stats.table";
      std::ifstream f(file_name.c_str());
      if (f) {
         std::string line;
         unsigned int n_grid_points = 0;
         std::map<unsigned int, std::tuple<double, double, double> > grid_map;
         while (std::getline(f, line)) {
            std::vector<std::string> words = coot::util::split_string_no_blanks(line);
            if (words.size() == 4) { // 5 with kurtosis
               unsigned int idx = util::string_to_int(words[0]);
               double mean = util::string_to_double(words[1]);
               double var  = util::string_to_double(words[2]);
               double skew = util::string_to_double(words[3]);
               std::tuple<double, double, double> t(mean, var, skew);
               grid_map[idx] = t;
            }
         }

         for (int ix= -n_steps; ix<=n_steps; ix++) {
            for (int iy= -n_steps; iy<=n_steps; iy++) {
               for (int iz= -n_steps; iz<=n_steps; iz++) {
                  int idx =
                     (ix + n_steps) * n_per_side * n_per_side +
                     (iy + n_steps) * n_per_side +
                     (iz + n_steps);
                  if (useable_grid_points.find(idx) != useable_grid_points.end()) {
                     clipper::Coord_orth pt_in_grid = make_pt_in_grid(ix, iy, iz, step_size, axes);
                     clipper::Coord_orth pt_grid_point = cb_pt + pt_in_grid;

                     std::cout << "check-stats mean "
                               << idx << " "
                               << pt_grid_point.x() << " "
                               << pt_grid_point.y() << " "
                               << pt_grid_point.z() << " "
                               << std::get<0>(grid_map[idx])
                               << std::endl;
                  }
               }
            }
         }
      }
   }
}


// static
void
coot::side_chain_densities::combine_directory(const std::string &rot_dir, int n_steps,
                                              double mn_unreliable_minimum_counts,
                                              double mn_unreliable_minimum_counts_for_low_variance,
                                              double mn_unreliable_minimum_variance,
                                              double mn_use_this_variance_for_unreliable) {

   // work out the mean and standard deviation of the grid point
   // for this rotamer (of this residue)

   // first check that the files are all standard format

   unsigned int n_per_side = 2 * n_steps + 1;
   unsigned int n_target = n_per_side * n_per_side * n_per_side;

   std::string glob_pattern = "*.tab";
   std::vector<std::string> files = coot::util::glob_files(rot_dir, glob_pattern);

   // first check that every file has the same number of numbers
   std::map<std::string, unsigned int> number_count_map;
   for (std::size_t i=0; i<files.size(); i++) {
      const std::string &fn = files[i];
      std::ifstream f(fn.c_str());
      if (f) {
         unsigned int word_count = 0;
         std::string line;
         while (std::getline(f, line)) {
            std::vector<std::string> words = coot::util::split_string_no_blanks(line);
            word_count += words.size();
         }
         number_count_map[fn] = word_count;
      } else {
         std::cout << "failed to open " << fn << std::endl;
      }
   }

   std::map<std::string, unsigned int>::const_iterator it;
   for (it=number_count_map.begin(); it!=number_count_map.end(); it++) {
      // std::cout << "counts " << it->first << " " << it->second << std::endl;
      if (it->second != n_target) {
         std::cout << "combine_directory() fail " << rot_dir << " " << it->first
                   << " " << it->second << " c.f. " << n_target << std::endl;
         return;
      }
   }


   // OK, files were good

   // std::cout << "debug:: files were good for " << rot_dir << std::endl;

   // so that we can organize the memory a bit easier - on the stack
   const int nps = n_steps * 2 + 1;
   const int nnn = nps * nps * nps;
   std::vector<float> *x = new std::vector<float>[nnn]; // deleted
   for (int i=0; i<nnn; i++) x[i].reserve(40);

   for (std::size_t i=0; i<files.size(); i++) {
      const std::string &fn = files[i];
      std::ifstream f(fn.c_str());
      if (f) {
         unsigned int word_count = 0;
         std::string line;
         while (std::getline(f, line)) {
            std::vector<std::string> words = coot::util::split_string_no_blanks(line);
            for (std::size_t ii=0; ii<words.size(); ii++) {
               const std::string &w = words[ii];
               try {
                  float f = util::string_to_float(w);
                  if (f > -1000.0) { // mask value
                     x[word_count].push_back(f);
                  }
               }
               catch (const std::runtime_error &rte) {
                  std::cout << "ERROR:: rte: combine_directory() " << fn << " " << rte.what() << std::endl;
               }
               word_count++;
            }
         }
         // is fine
         // std::cout << "ending word count "<< word_count << " " << fn << std::endl;
      }
   }

   if (true) {

      stats::single s_all;

      for (int i=0; i<nnn; i++) {

         // stats
         {
            if (x[i].size() > 1) {
               stats::single s;
               for (std::size_t j=0; j<x[i].size(); j++) {
                  const float &val = x[i][j];
                  s.add(val);
                  s_all.add(val);
               }

               std::string fn_stats =  + "stats.table";
               fn_stats = util::append_dir_file(rot_dir, fn_stats);
               std::ios_base::openmode mode = std::ios_base::app;
               std::ofstream f_stats(fn_stats.c_str(), mode);
               double mean = s.mean();
               double var  = s.variance();
               double skew = s.skew();
               std::pair<double, double> mi = s.median_and_iqr();
               double median = mi.first;
               double irq    = mi.second;
               double lim_low  = median - 1.75 * irq;
               double lim_high = median + 1.75 * irq;

               // remove values that are more than 3 iqrs from the mean
               stats::single s_filtered;
               for (std::size_t j=0; j<x[i].size(); j++) {
                  const float &val = x[i][j];
                  if (val < lim_high)
                     if (val > lim_low)
                        s_filtered.add(val);
               }

               bool unreliable = false;
               unsigned int n_filtered = s.size() - s_filtered.size();
               std::cout << "debug:: n_filtered " << n_filtered << std::endl;
               if (s_filtered.size() < mn_unreliable_minimum_counts) // magic number 
                  unreliable = true;
               if (s_filtered.size() < mn_unreliable_minimum_counts_for_low_variance &&  // magic numbers
                   var < mn_unreliable_minimum_counts_for_low_variance)
                  unreliable = true;

               if (unreliable)
                  var = mn_use_this_variance_for_unreliable; // magic number

               if (f_stats) {
                  f_stats << i << " " << s_filtered.mean() << " " << s_filtered.variance()
                          << " " << s_filtered.skew() << "\n";
               } else {
                  std::cout << "failed to open " << fn_stats << std::endl;
               }
               f_stats.close();

            } else {

               if (x[i].size() == 1) {
                  std::cout << "only 1 point " << rot_dir << " " << i << std::endl;
                  double mean = x[i][0];
                  double var  = 1.0;
                  double skew = 0.0;
                  std::string fn_stats =  + "stats.table";
                  fn_stats = util::append_dir_file(rot_dir, fn_stats);
                  std::ios_base::openmode mode = std::ios_base::app;
                  std::ofstream f_stats(fn_stats.c_str(), mode);
                  f_stats << i << " " << mean << " " << var << " " << skew << "\n";

               } else {
                  // std::cout << "No grid point data for " << rot_dir << " grid point " << i << std::endl;
               }
            }
         }

         // file output
         bool file_output = true; // maybe not... (useful for R, not otherwise and makes many files)
         file_output = false; // Not at the moment, then.
         if (file_output) {
            std::string fn = "grid-point-" + util::int_to_string(i) + ".data";
            fn = util::append_dir_file(rot_dir, fn);
            if (x[i].size() > 0) {
               std::ofstream f(fn.c_str());
               if (f) {
                  for (std::size_t j=0; j<x[i].size(); j++)
                     f << x[i][j] << "\n";
                  f.close();
               } else {
                  std::cout << "Failed to open " << fn << " for writing " << std::endl;
               }
            } else {
               // std::cout << "No grid point data for " << rot_dir << " grid point " << i << std::endl;
            }
         }

         bool screen_output = false;
         if (screen_output) {
            std::cout << "sample point of nnn: " << i << " " << x[i].size() << " values" << std::endl;
            for (std::size_t j=0; j<x[i].size(); j++) {
               std::cout << x[i][j] << " ";
            }
            std::cout << std::endl;
         }
      }

      // all grid points, all sample data:
      std::string fn_summary_stats = "summary-stats";
      fn_summary_stats = util::append_dir_file(rot_dir, fn_summary_stats);
      // std::ios_base::openmode mode = std::ios_base::fwrite;
      std::ofstream f_stats(fn_summary_stats.c_str());
      double mean_all = s_all.mean();
      double var_all = s_all.variance();

      f_stats << s_all.size() << " " << mean_all << " " << var_all << std::endl;
      f_stats.close();

   }

   delete [] x;

}


void
coot::side_chain_densities::set_magic_number(const std::string &mn_name, double val) {

   get_results_addition_lock();
   if (mn_name == "mn_log_likelihood_ratio_difference_min") mn_log_likelihood_ratio_difference_min = val;
   if (mn_name == "mn_scale_for_normalized_density") mn_scale_for_normalized_density = val;
   if (mn_name == "mn_density_block_sample_x_max") mn_density_block_sample_x_max = val;
   release_results_addition_lock();

}
