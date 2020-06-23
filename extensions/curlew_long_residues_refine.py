

def long_residues_refine(n_residues=3):
    active_atom = closest_atom_simple_py() # active_atom returns the CA if it can
    if not active_atom:
       print "No active atom"
    else:
       imol       = active_atom[0]
       chain_id   = active_atom[1]
       res_no     = active_atom[2]
       ins_code   = active_atom[3]
       atom_name  = active_atom[4]
       alt_conf   = active_atom[5]
       specs = []
       for ires in range(res_no-n_residues, res_no+n_residues+1):
           # test if the residue exists by looking for a residue name
           rn = residue_name(imol, chain_id, ires, ins_code)
           if len(rn) > 0:
               specs.append([chain_id, ires, ins_code])
       refine_residues(imol, specs)



register_extension("curlew_long_residues_refine.py", "1.0")
