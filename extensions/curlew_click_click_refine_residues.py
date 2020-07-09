
def click_click_refine_residues():

    def click_click_func(*clicks):
        if len(clicks) == 2:
            # has to be!?
            click_1 = clicks[0]
            click_2 = clicks[1]
            print "click_1:", click_1
            print "click_2:", click_2
            if ((len(click_1) == 7) and
                (len(click_2) == 7)):
                resname_1 = residue_name(*click_1[1:5])

                imol_click_1 = click_1[1]
                imol_click_2 = click_2[1]

                chain_id_1 = click_1[2]
                chain_id_2 = click_2[2]

                res_no_1 = click_1[3]
                res_no_2 = click_2[3]

                if not imol_click_1 == imol_click_2:
                    # some error message?
                    return False
                else:
                    imol = imol_click_1
                    chain_id = chain_id_1

                    if not chain_id_1 == chain_id_2:
                        return False
                    else:
                        if res_no_1 > res_no_2:
                            res_no_1, res_no_2 = res_no_2, res_no_1
                        numbers = range(res_no_1, res_no_2 + 1)
                        residue_specs = map(lambda rn: [chain_id, rn, ""],
                                            numbers)
                        refine_residues(imol, residue_specs)

    user_defined_click(2, click_click_func)

if coot_python.main_menubar:

    jligand_menu = coot_menubar_menu("Refine")

    add_simple_coot_menu_menuitem(
        jligand_menu, "Click-click Refine Residues",
        lambda func: click_click_refine_residues())

register_extension("curlew_click_click_refine_residues.py", "1.1")
