
def chain_refine(imol, chain_id):
    if not valid_map_molecule_qm(imol_refinement_map()):

        info_dialog("WARNING:: Must set the refinement map")

        ls = residues_in_chain(imol, chain_id)
        refine_residues(imol, ls)

def black_box_morph_and_fit(imol, chain_id, imol_map):

    save_matrix = matrix_state()

    set_matrix(3 * save_matrix)
    set_use_trans_peptide_restraints(0)

    generate_local_self_restraints(imol, chain_id, 4.1)
    set_refinement_geman_mcclure_alpha(1.0)
    set_show_extra_restraints(imol, 0)

    with AutoAccept():
        morph_fit_by_secondary_structure_elements(imol, chain_id)
        chain_refine(imol, chain_id)

    set_matrix(save_matrix)
    set_refinement_geman_mcclure_alpha(0.05)

    with AutoAccept():
        chain_refine(imol, chain_id)

    set_refinement_geman_mcclure_alpha(0.2)

def ammof(imol, chain_id, imol_map, pir_alignment_file):
    spf = dragged_refinement-steps_per_frame()
    set_dragged_refinement_steps_per_frame(300)
    associate_pir_alignment_from_file(imol, chain_id, pir_alignment_file)
    apply_pir_alignment(imol, chain_id)
    black_box_morph_and_fit(imol, chain_id, imol_map)
    set_dragged_refinement_steps_per_frame(spf)

register_extension("curlew_black_box_morph_and_fit.py", "1.1")
