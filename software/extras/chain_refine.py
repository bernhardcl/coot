
def residues_in_chain(imol, chain_id_in):
    res = residues_matching_criteria(imol, lambda chain_id, res_no, ins_code,
                                     serial_no: (chain_id == chain_id_in))
    return res

def chain_ref_func(refine=True):

    if not valid_map_molecule_qm(imol_refinement_map()):
        info_dialog("Must set the refinement map.")
    else:
        with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf]:

            rc_spec = [aa_chain_id, aa_res_no, aa_ins_code]
            ls = residues_in_chain(aa_imol, aa_chain_id)

            if refine:
                refine_residues(aa_imol, ls)
            else:
                # regularize
                regularize_residues(aa_imol, ls)

add_key_binding("Chain Refinement", "E", lambda: chain_ref_func())

add_key_binding("Chain Regularize", "Y", lambda: chain_ref_func(False))
