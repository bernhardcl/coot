
try:
    # to see if it is defined
    pepflip_intermediate_atoms
    add_key_binding("Pepflip Complicado", "q",
                    lambda:
                    pepflip_intermediate_atoms() == 0 and
                    pepflip_active_residue())
except:
    pass

register_extension("curlew_pepflip_complicado.py", 1.1)
