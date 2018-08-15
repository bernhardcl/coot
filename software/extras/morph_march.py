if (have_coot_python):
    if coot_python.main_menubar():
        menu = coot_menubar_menu("Morph")

        add_simple_coot_menu_menuitem(
            menu,
            "Morph Fit This Molecule",
            lambda func:
            using_active_atom(morph_fit_all, "aa_imol", 10)
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Morph Fit This Chain 8",
            lambda func:
            using_active_atom(morph_fit_chain, "aa_imol", "aa_chain_id", 8)
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Morph Fit This Chain 12",
            lambda func:
            using_active_atom(morph_fit_chain, "aa_imol", "aa_chain_id", 12)
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Morph Fit This Chain 15",
            lambda func:
            using_active_atom(morph_fit_chain, "aa_imol", "aa_chain_id", 15)
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Morph Fit This Chain 18",
            lambda func:
            using_active_atom(morph_fit_chain, "aa_imol", "aa_chain_id", 18)
            )


        # def morph_fit_region_func(select_radius, morph_radius):
        #     with UsingActiveAtom() as \
        #     [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
        #      aa_atom_name, aa_alt_conf]:
        #         central_residue = [aa_chain_id, aa_res_no, aa_ins_code]
        #         other_residues = residues_near_residue(aa_imol,
        #                                                central_residue,
        #                                                select_radius)
        #         morph_residues = [central_residue] + other_residues
        #         print "BL DEBUG:: morph_res", morph_residues
        #         morph_fit_residues(aa_imol, morph_residues, morph_radius)

        # add_simple_coot_menu_menuitem(
        #     menu,
        #     "Morph Fit This Region 12 9",
        #     lambda func:
        #     morph_fit_region_func(12, 9)
        #     )

        # add_simple_coot_menu_menuitem(
        #     menu,
        #     "Morph Fit This Region 8 6",
        #     lambda func:
        #     morph_fit_region_func(8, 6)
        #     )


        add_simple_coot_menu_menuitem(
            menu,
            "Secondary Structure Morph Fit This Chain",
            lambda func:
            using_active_atom(morph_fit_by_secondary_structure_elements,
                              "aa_imol", "aa_chain_id")
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Jiggle-fit This Molecule",
            lambda func:
            using_active_atom(fit_molecule_to_map_by_random_jiggle,
                              "aa_imol", 1000, 2)
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Jiggle-fit This Chain",
            lambda func:
            using_active_atom(fit_chain_to_map_by_random_jiggle,
                              "aa_imol", "aa_chain_id", 2000, 2)
            )


        def jiggle_all_chains_func():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                for chain_id in chain_ids(aa_imol):
                    fit_chain_to_map_by_random_jiggle(aa_imol, chain_id, 800, 2)
                    fit_chain_to_map_by_random_jiggle(aa_imol, chain_id, 800, 2)

        add_simple_coot_menu_menuitem(
            menu,
            "Jiggle All Chains",
            lambda func: jiggle_all_chains_func()
            )

