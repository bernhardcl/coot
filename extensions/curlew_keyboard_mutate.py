
#
# type is a single letter code string
def keyboard_mutate_active_residue(residue_type):

    print "mutate to this residue_type:", residue_type

    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:

        dic_type_1lc = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
        # NOTE: probably should use a dictionary the other way around!?
        for tlc, slc in dic_type_1lc.items():
            if residue_type.upper() == slc:
                print "mutate this: aa-res-no: %s to %s" %(aa_res_no, tlc)
                mutate_and_auto_fit(aa_res_no, aa_chain_id, aa_imol,
                                    imol_refinement_map(), tlc)

def mutate_and_auto_fit_window():

    def delete_event(*args):
        window.destroy()
        return False

    def do_entry_event(widget, event):
        if (event.keyval == 65293):  # return
            text = entry.get_text()
            keyboard_mutate_active_residue(text)
            window.destroy()
        return False

    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    vbox = gtk.VBox(False, 2)
    message_text = gtk.Label("Mutate and Autofit")
    h_sep = gtk.HSeparator()
    hbox_for_buttons = gtk.HBox(False, 4)
    entry = gtk.Entry()
    close_button = gtk.Button("   Close   ")

    window.set_title("Mutate and Autofit")
    vbox.pack_start(message_text, False, False, 4)
    vbox.pack_start(h_sep, False, False, 4)
    vbox.pack_start(hbox_for_buttons, False, False, 4)
    hbox_for_buttons.pack_start(entry, False, False, 4)
    vbox.pack_start(close_button, False, False, 4)
    hbox_for_buttons.set_homogeneous(True)
    window.add(vbox)

    entry.connect("key_press_event", do_entry_event)
    close_button.connect("clicked", delete_event)
    window.show_all()

add_key_binding("Fast auto-mutate", "q",
                lambda: mutate_and_auto_fit_window())

register_extension("curlew_keyboard_mutate.py", "1.0")
