
"""
@author: nabin

- Validates labeling of Cryo2Struct dataset. 

"""

import os
import mrcfile
from copy import deepcopy
from Bio import PDB

import warnings
warnings.filterwarnings("ignore")


def get_xyz(idx, voxel, origin):
    return (idx * voxel) + origin


def check_elements_in_list(idx_xyzs, known_struct_coordinates):

    coordinate_not_found = list()
    for element in idx_xyzs:
        if element not in known_struct_coordinates:
            coordinate_not_found.append(coordinate_not_found)

    return coordinate_not_found


def extract_known_coordinates(known_structure):
    parser = PDB.PDBParser()
    structure = parser.get_structure('protein', known_structure)

    known_structure_coordinate_list = list()

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "CA":
                        x, y, z = atom.get_coord()
                        known_structure_coordinate_list.append([int(x), int(y), int(z)])
    return known_structure_coordinate_list


def idx_xyz_conversion(mrc_map):
    x_origin = mrc_map.header.origin['x']
    y_origin = mrc_map.header.origin['y']
    z_origin = mrc_map.header.origin['z']
    x_voxel = mrc_map.voxel_size['x']
    y_voxel = mrc_map.voxel_size['y']
    z_voxel = mrc_map.voxel_size['z']
    idx_xyz_list = list()
    mrc_data = deepcopy(mrc_map.data)
    index_error = False
    for k in range(len(mrc_data[2])):
        for j in range(len(mrc_data[1])):
            for i in range(len(mrc_data[0])):
                try:
                    if mrc_data[i][j][k] > 0:
                        x = round(get_xyz(k, x_voxel, x_origin), 3)
                        y = round(get_xyz(j, y_voxel, y_origin), 3)
                        z = round(get_xyz(i, z_voxel, z_origin), 3)
                        idx_xyz_list.append([int(x), int(y), int(z)])
                except IndexError:
                    index_error = True
    print("Index Error :", index_error)

    return idx_xyz_list



if __name__ == "__main__":
    input_directory = "cryo2struct/root/EMD_1"
    density_maps = os.listdir(input_directory)
    for maps in density_maps:
        maps_dir = f"{input_directory}/{maps}"
        input_map = f"{maps_dir}/atom_ca_emd_normalized_map.mrc"
        pdb_files = [f for f in os.listdir(maps_dir) if f.endswith(".pdb")][0]
        pdb_file_name = pdb_files.split(".")[0]
        pdb_file_name = pdb_file_name.split("_")[0]
        known_pdb_backbone = f"{maps_dir}/{pdb_file_name}.pdb"

        mrc_map = mrcfile.open(input_map, mode='r')
        idx_xyzs=idx_xyz_conversion(mrc_map)
        known_struct_coordinates = extract_known_coordinates(known_pdb_backbone)
        print(f"Residues == Map : {maps}, Conversion : {len(idx_xyzs)},  Known Structure: {len(known_struct_coordinates)}")
        cord_not_found = check_elements_in_list(idx_xyzs, known_struct_coordinates)   
        print(f"Coordinates of conversion not found in known structure of {maps} : {len(cord_not_found)}") 
        print("######################################################################################")





