def rama_outlier_gui():
    """
    A gui to list Ramachandran outliers etc.
    A first draft, may become more sophisticated at some point
    """

    def list_rama_outliers(imol):

        r = all_molecule_ramachandran_region(imol)
        outliers = []
        allowed = []
        for res in r:
            if res[1] == 0:
                outliers.append(res[0])
            if res[1] == 1:
                allowed.append(res[0])

        def make_buttons(res_list, label_string):
            ret = []
            for res_spec in res_list:
                chain_id = res_spec[1]
                res_no = res_spec[2]
                label = label_string + ": " + \
                        chain_id + " " + str(res_no)
                func = [cmd2str(set_go_to_atom_molecule, imol),
                        cmd2str(set_go_to_atom_from_res_spec, res_spec)]
                ret.append([label, func])
                
            return ret

        outlier_buttons = make_buttons(outliers, "Outlier")
        allowed_buttons = make_buttons(allowed, "Allowed")
        all_buttons = outlier_buttons + allowed_buttons
        
        def clear_and_add_back(vbox, outliers_list, allowed_list, filter_flag):
            # clear
            children = vbox.get_children()
            map(lambda c: c.destroy(), children)
            # add back
            if not filter_flag:
                buttons = outliers_list + allowed_list
            else:
                # filter
                buttons = outliers_list
            map(lambda button_info: add_button_info_to_box_of_buttons_vbox(button_info, vbox),
                buttons)

        dialog_box_of_buttons_with_check_button(
            " Ramachandran issues ", [300, 300], [], "  Close  ",
            "Outliers only",
            lambda check_button, vbox: clear_and_add_back(vbox, outlier_buttons, allowed_buttons, True)
            if check_button.get_active() else
            clear_and_add_back(vbox, outlier_buttons, allowed_buttons, False),
            False)

   
    molecule_chooser_gui(
      "List Rama outliers for which molecule?",
          lambda imol: list_rama_outliers(imol))

menu = coot_menubar_menu("Validate")
if menu:
    add_simple_coot_menu_menuitem(menu, "List Ramachandran outliers...",
                                  lambda func: rama_outlier_gui())

# rama_outlier_gui()
