
try:
    pepflip_intermediate_atoms
    add_key_binding("Pepflip Complicado", "q",
                    lambda func:
                    pepflip_active_residue() if pepflip_intermediate_atoms() == 0)
except:
    pass

register_extension("curlew_pepflip_complicado.py", 1.1)
