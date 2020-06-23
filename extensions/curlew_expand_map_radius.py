
def curlew_map_radius_change(increase_flag):
    ml = map_molecule_list()
    imol_map = ml[0]
    if True:  # at some point: is_EM_map_qm(imol_map)
        r = get_map_radius()
        set_map_radius(r * 1.5 if increase_flag else r * 0.6666665)

add_key_binding("Expand map radius", "]", lambda: curlew_map_radius_change(True))

add_key_binding("Decrease map radius", "[", lambda: curlew_map_radius_change(False))

register_extension("curlew_expand_map_radius.py", "1.0")
