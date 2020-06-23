# let's try to add neighbour-refine and next residue recentering
# 
def key_binding_terminal_residue_complicado():
    active_atom = active_residue()
    if (not active_atom):
        print "No active atom"
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]

        # which terminus is this?
        terminus_type='U' # unknown
        atoms_prev = residue_info(imol, chain_id, res_no-1, ins_code)
        atoms_next = residue_info(imol, chain_id, res_no+1, ins_code)

        # print "atoms_prev", atoms_prev
        # print "atoms_next", atoms_next

        if atoms_next == False:
            terminus_type = 'C'
            new_res_no = res_no + 1
            
        if atoms_prev == False:
            terminus_type = 'N'
            new_res_no = res_no - 1

        set_dragged_refinement_steps_per_frame(520)
        set_terminal_residue_do_rigid_body_refine(0)
        set_add_terminal_residue_n_phi_psi_trials(10000); # default 1000
        # set_add_terminal_residue_add_other_residue_flag(1)
        success = add_terminal_residue(imol, chain_id, res_no, "auto", 1)

        if success == 1:
            if terminus_type != 'U':
                
                rc_spec = [chain_id, res_no, ins_code]
                ls = residues_near_residue(imol, rc_spec, 1.9)
                with_auto_accept([refine_residues, imol, [rc_spec] + ls])

                set_go_to_atom_molecule(imol)
                set_go_to_atom_chain_residue_atom_name(chain_id, new_res_no, ' CA ');

def key_binding_delete_N_terminal_of_active_residue():
    active_atom = active_residue()
    if not active_atom:
       print("No active atom")
    else:
       imol = active_atom[0]
       residue_spec = atom_spec_to_residue_spec(active_atom)
       del_residue_spec = residue_spec # initially
       N_terminal_residue_spec = [residue_spec_to_chain_id(residue_spec),
                                  residue_spec_to_res_no(residue_spec)-1,
                                  residue_spec_to_ins_code(residue_spec)]
       # not to self: residue_info should work with a spec
       ri = residue_info(imol,
                         residue_spec_to_chain_id(N_terminal_residue_spec),
                         residue_spec_to_res_no(N_terminal_residue_spec),
                         residue_spec_to_ins_code(N_terminal_residue_spec))
       try:
          l = len(ri)
          del_residue_spec = N_terminal_residue_spec
          delete_residues(imol, [del_residue_spec])
       except:
          # no next residue
          delete_residues(imol, [del_residue_spec])

def key_binding_refine_triple():          
    active_atom = active_residue()
    if not active_atom:
       print("No active atom")
    else:
       imol = active_atom[0]
       residue_spec = atom_spec_to_residue_spec(active_atom)
       N_terminal_residue_spec = [residue_spec_to_chain_id(residue_spec),
                                  residue_spec_to_res_no(residue_spec)-1,
                                  residue_spec_to_ins_code(residue_spec)]
       C_terminal_residue_spec = [residue_spec_to_chain_id(residue_spec),
                                  residue_spec_to_res_no(residue_spec)+1,
                                  residue_spec_to_ins_code(residue_spec)]
       spec_list = [N_terminal_residue_spec, residue_spec, C_terminal_residue_spec]
       refine_residues(imol, spec_list)

def key_binding_terminal_spin():
    active_atom = active_residue()
    if not active_atom:
       print("No active atom")
    else:
       imol = active_atom[0]
       residue_spec = atom_spec_to_residue_spec(active_atom)
       print ('spin_N {} {} {}'.format(imol, residue_spec, 120))
       spin_N_py(imol, residue_spec, 120)


# add_key_binding("Add terminal residue", "y", lambda: key_binding_func_1())
add_key_binding("Add terminal residue", "y", lambda: key_binding_terminal_residue_complicado())
add_key_binding("Spin N-and-sidechain", "Y", lambda: key_binding_terminal_spin())
add_key_binding("Delete Ative Resiude", "!", lambda: key_binding_delete_N_terminal_of_active_residue())
add_key_binding("Refine triple", "@", lambda: key_binding_refine_triple())

register_extension("curlew_add_terminal_residue_complicado.py", "1.0")
