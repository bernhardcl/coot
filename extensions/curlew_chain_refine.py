
try:
    # check if we have lined_residues function
    linked_residues_py
    def refine_active_fragment():
        with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf,
                                       aa_res_spec]:
            residues = linked_residues_py(aa_res_spec, aa_imol, 1.7)
            if isinstance(residues, list):
                refine_residues(aa_imol, residues)
            else:
                print "residues is not a list!"
except:
    refine_active_fragment = False

def add_self_restraints_for_this_glyco_tree():

    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        gtr = glyco_tree_residues(aa_imol, aa_res_spec)
        generate_local_self_restraints_by_residues(aa_imol, gtr, 5.5)

def refine_this_glyco_tree():

    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        gtr = glyco_tree_residues(aa_imol, aa_res_spec)
        refine_residues(aa_imol, gtr)


def chain_refine():
    if not valid_map_molecule_qm(imol_refinement_map()):
        info_dialog("Must set the refinement map")
    else:
        with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf]:
            ls = residues_in_chain(aa_imol, aa_chain_id)
            refine_residues(aa_imol, ls)

add_key_binding("Chain Refinement", "E", lambda: chain_refine())


if (have_coot_python):
  if coot_python.main_menubar():
      menu = coot_menubar_menu("Refine")

      def refine_all_atoms_func():
          with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                     aa_ins_code, aa_atom_name, aa_alt_conf]:
              ls = all_residues(aa_imol)
              refine_residues(aa_imol, ls)

      add_simple_coot_menu_menuitem(
          menu, "All-atom Refine",
          lambda func: refine_all_atoms_func())

      add_simple_coot_menu_menuitem(
          menu, "Chain Refine",
          lambda func: chain_refine())

      if refine_active_fragment:
          add_simple_coot_menu_menuitem(
              menu, "Refine Fragment",
              lambda func: refine_active_fragment())

      add_simple_coot_menu_menuitem(
          menu, "Make Glyco Tree Self Restraints",
          lambda func: add_self_restraints_for_this_glyco_tree())

      add_simple_coot_menu_menuitem(
          menu, "Refine This Glyco Tree",
          lambda func: refine_this_glyco_tree())

      for w in [0.01, 0.03, 0.1, 0.3, 1, 3, 10]:
          add_simple_coot_menu_menuitem(
              menu, "Set Geman-McClure Alpha to " + str(w),
              lambda func: set_refinement_geman_mcclure_alpha(w))


register_extension("curlew_chain_refine.py", "1.3")
