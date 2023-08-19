
"""
@author: nabin

- Validates the density maps, MRC format.
"""


import mrcfile
import os


def validate_mrc_file(mrc_file_path):
    density_map_mrc = mrcfile.open(mrc_file_path, mode='r+')
    density_map_mrc.update_header_stats()
    density_map_mrc.close()
    try:
        val = mrcfile.validate(name=mrc_file_path)
        if val == False:
            print(f"File does not meet the MRC format specification in any way \n {mrc_file_path}")
            exit()
    except RuntimeWarning:
        print(f"File is seriously invalid \n {mrc_file_path}")
        exit()

if __name__ == "__main__":
    input_directory = "cryo2struct/root/EMD_0"
    density_maps = os.listdir(input_directory)
    for maps in density_maps:
        normalized_map = f"{input_directory}/{maps}/emd_normalized_map.mrc"
        validate_mrc_file(normalized_map)
        ca_atom_label_map = f"{input_directory}/{maps}/atom_ca_emd_normalized_map.mrc"
        validate_mrc_file(ca_atom_label_map)
        atom_label_map = f"{input_directory}/{maps}/atom_emd_normalized_map.mrc"
        validate_mrc_file(atom_label_map)
        amino_label_map = f"{input_directory}/{maps}/amino_emd_normalized_map.mrc"
        validate_mrc_file(amino_label_map)
        sec_label_map = f"{input_directory}/{maps}/sec_struc_emd_normalized_map.mrc"
        validate_mrc_file(sec_label_map)
    print("ALL DONE!!", input_directory)
