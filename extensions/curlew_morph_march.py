
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
                              "aa_imol", 1500, 2)
            )


        add_simple_coot_menu_menuitem(
            menu,
            "Jiggle-fit This Chain",
            lambda func:
            using_active_atom(fit_chain_to_map_by_random_jiggle,
                              "aa_imol", "aa_chain_id", 3000, 2.2)
            )


        def x():
            pass

        add_simple_coot_menu_menuitem(
            menu,
            "Jiggle All Chains",
            lambda func: map(using_active_atom(fit_chain_to_map_by_random_jiggle, "aa_imol", chain_id), chain_ids("aa_imol"))
            )



