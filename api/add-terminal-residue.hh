#ifndef ADD_TERMINAL_RESIDUE_HH
#define ADD_TERMINAL_RESIDUE_HH

#include <string>
#include <utility>
#include "compat/coot-sysdep.h"
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include "coot-utils/atom-selection-container.hh"
#include "utils/ctpl.h"

namespace coot {

   void remove_TER_internal(mmdb::Manager *mol, mmdb::Residue *residue_p);
   void remove_OXT_internal(mmdb::Residue *residue_p, mmdb::Manager *mol); // usually not there, but if it is, it needs to go.

   atom_selection_container_t
   add_side_chain_to_terminal_res(atom_selection_container_t asc,
                                  const std::string &res_type,
                                  const std::string &terminus_type,
                                  float b_factor_for_new_atoms,
                                  const protein_geometry &geom);

   std::pair<int, std::string>
   add_terminal_residue(int imol_no, const std::string &terminus_type, mmdb::Residue *residue_p,
                        mmdb::Manager *mol, int udd_atom_index,
                        const std::string &chain_id, const std::string &res_type, float b_factor_for_new_atoms,
                        const clipper::Xmap<float> &xmap, const protein_geometry &geom,
                        ctpl::thread_pool &static_thread_pool);

   void remove_ter_atoms(mmdb::Manager *mol, mmdb::Residue *residue_p);

   std::string get_term_type(mmdb::Residue *residue_p, mmdb::Manager *mol);

}

#endif // ADD_TERMINAL_RESIDUE_HH
