
ligand_check_refmac_mtz = False
ligand_check_refmac_fobs_col = False
ligand_check_refmac_sigfobs_col = False
ligand_check_refmac_r_free_col = False

def ligand_check_refmac_columns(f_list, sigf_list, rfree_list):

    # dummy function, not needed anywhere (as is)
    # need the above as globals...

    # happy path

    # Using the first sigf (there should only be one typically)
    # get the F label (Fx) from x/y/SIGFPx
    #
    pass

def jiggle_fit_active_residue():
    with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        coot.fit_to_map_by_random_jiggle(aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, 100, 1.0)

if (use_gui_qm != 2):
    menu = coot_gui.coot_menubar_menu("Ligand")

    coot_gui.add_simple_coot_menu_menuitem(
        menu,
        "Find Ligands...",
        lambda func: coot.do_find_ligands_dialog())

    coot_gui.add_simple_coot_menu_menuitem(
        menu,
        "Jiggle-Fit Ligand",
        lambda func: jiggle_fit_active_residue())

    coot_gui.add_simple_coot_menu_menuitem(
        menu,
        "Hydrogenate region",
        lambda func: coot.hydrogenate_region(6))

    coot_gui.add_simple_coot_menu_menuitem(
          menu,
          "SMILES -> 2D",
          lambda func:
          coot_gui.generic_single_entry("SMILES string",
                               "", " Send to 2D Viewer ",
                               lambda text: coot.smiles_to_ligand_builder(text)))

    coot_gui.add_simple_coot_menu_menuitem(
          menu,
          "SMILES -> simple 3D",
          lambda func:
          coot_gui.generic_double_entry("Residue name", "SMILES string  ", "LIG", "",
                               False, False, 
                               "Import Molecule",
                               lambda text_1, text_2:
                               coot.import_rdkit_mol_from_smiles(text_1, text_2)))


    coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Residue -> 2D",
            lambda func:
            coot_utils.using_active_atom(residue_to_ligand_builder,
                              "aa_imol", "aa_chain_id", "aa_res_no",
                              "aa_ins_code", 0.015))

    if coot.enhanced_ligand_coot_p():
        def flev_rdkit_func():
            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                prodrg_import.fle_view_with_rdkit(aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, 4.2)
                coot.set_flev_idle_ligand_interactions(1)
            
        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "FLEV this residue",
            lambda func: flev_rdkit_func()
            )


        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Toggle FLEV Ligand Interactions",
            lambda func: toggle_idle_ligand_interactions()
            )

        def go_solid_func(state):
            set_display_generic_objects_as_solid(state)
            graphics_draw()

        add_simple_coot_menu_menuitem(
             menu, "Solid Generic Objects",
             lambda func: go_solid_func(1))

        add_simple_coot_menu_menuitem(
             menu, "Unsolid Generic Objects",
             lambda func: go_solid_func(0))


        def show_chem_func():
            set_display_generic_objects_as_solid(1) # there may be consequences...
            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                coot.show_feats(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)
            
        coot_gui.add_simple_coot_menu_menuitem(
            menu,
            "Show Chemical Features",
            lambda func: show_chem_func()
            )

    def rename_atoms_to_reference(menuitem):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                   aa_ins_code, aa_atom_name, aa_alt_conf]:
           # user chooses cif_dict_file_name_out (provide default)
           # reference_comp_id and new_comp_id

           # args: entry_info_list, check_button_info, go_button_label, handle_go_function
           coot_gui.generic_multiple_entries_with_check_button(
               [["      Reference Residue Type: ", "ATP"],
                ["               New Residue Type: ", "RXC"],
                ["New Dictionary cif file name: ", "new-dictionary-RXC.cif"]],
               False,
               "   OK   ",
               lambda ref_comp_id_txt_new_comp_id_txt_new_cif_file_name_txt: coot.match_this_residue_and_dictionary(aa_imol, aa_chain_id, aa_res_no, aa_ins_code, ref_comp_id_txt_new_comp_id_txt_new_cif_file_name_txt[2], ref_comp_id_txt_new_comp_id_txt_new_cif_file_name_txt[0], ref_comp_id_txt_new_comp_id_txt_new_cif_file_name_txt[1]))
        

    coot_gui.add_simple_coot_menu_menuitem(
        menu, "Rename Atom to Reference", rename_atoms_to_reference)

    def go_solid_func(state):
        set_display_generic_objects_as_solid(state)
        coot.graphics_draw()

    coot_gui.add_simple_coot_menu_menuitem(
         menu, "Solid Generic Objects",
         lambda func: go_solid_func(1))

    coot_gui.add_simple_coot_menu_menuitem(
         menu, "Unsolid Generic Objects",
         lambda func: go_solid_func(0))


    def tab_ligand_distortions_func():
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                   aa_ins_code, aa_atom_name, aa_alt_conf]:
             coot.print_residue_distortions(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)

    def display_ligand_distortions_func():
         with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
             set_display_generic_objects_as_solid(1)
             coot.display_residue_distortions(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)

    coot_gui.add_simple_coot_menu_menuitem(
         menu,
         "Tabulate (on terminal) Ligand Distorsions",
         lambda func: tab_ligand_distortions_func()
         )

    coot_gui.add_simple_coot_menu_menuitem(
         menu,
         "Display Ligand Distortions",
         lambda func: display_ligand_distortions_func())

    add_simple_coot_menu_menuitem(
        menu,
        "Quick Ligand Validate",
        lambda func: gui_ligand_check_dialog_active_residue())


##        # not interesting for the normal user!?
##        def density_ligand_score_func():
##            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
##                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
##               spec = [aa_chain_id, aa_res_no, aa_ins_code]
##               r = coot.density_score_residue(aa_imol, spec, coot.imol_refinement_map())
##               # BL says:: maybe a dialog for this?!
##               print "density at ligand atoms:", r

##        coot_gui.add_simple_coot_menu_menuitem(
##            menu,
##            "Density Score Ligand",
##            lambda func: density_ligand_score_func()
##            )


##        def fetch_ligand_pdbe_func():
##            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
##                               aa_ins_code, aa_atom_name, aa_alt_conf]:
##               comp_id = coot.residue_name(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)
##               print "here with residue name", comp_id
##               s = coot_utils.get_SMILES_for_comp_id_from_pdbe(comp_id)
##               if (isinstance(s, str)):
##                   pdbe_cif_file_name = os.path.join("coot-download", "PDBe-" + comp_id + ".cif")
##                   prodrg_import.import_from_3d_generator_from_mdl(pdbe_cif_file_name, comp_id)

##        coot_gui.add_simple_coot_menu_menuitem(
##            menu,
##            "### [Fetch ligand description & generate restraints]",
##            lambda func: fetch_ligand_pdbe_func()
##            )


    def probe_ligand_func():
        global probe_command
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                            aa_ins_code, aa_atom_name, aa_alt_conf]:
            ss = "//" + aa_chain_id + "/" + str(aa_res_no)
            imol_selection = coot.new_molecule_by_atom_selection(aa_imol, ss)
            work_dir = coot_utils.get_directory("coot-molprobity")
            tmp_selected_ligand_for_probe_pdb = os.path.join(work_dir,
                                                "tmp-selected-ligand-for-probe.pdb")
            tmp_protein_for_probe_pdb = os.path.join(work_dir,
                                                 "tmp-protein-for-probe.pdb")
            generic_objects.probe_dots_file_name = os.path.join(work_dir, "probe.dots")

            coot.set_mol_displayed(imol_selection, 0)
            coot.set_mol_active(imol_selection, 0)
            # BL comment: we assume H and view is correct.
            #set_go_to_atom_molecule(imol)
            #rc = residue_centre(imol, chain_id, res_no, ins_code)
            #set_rotation_centre(*rc)
            #hydrogenate_region(6)
            coot.write_pdb_file(imol_selection, tmp_selected_ligand_for_probe_pdb)
            coot.write_pdb_file(imol, tmp_protein_for_probe_pdb)
            coot_utils.popen_command(probe_command, ["-u", "-once", str(aa_res_no), # -once or -both
                                             "not " + str(aa_res_no),
                                             "-density60",
                                             tmp_selected_ligand_for_probe_pdb,
                                             tmp_protein_for_probe_pdb],
                                             [], generic_objects.probe_dots_file_name, False)
            coot.handle_read_draw_probe_dots_unformatted(dots_file_name, aa_imol, 0)
            coot.graphics_draw()



