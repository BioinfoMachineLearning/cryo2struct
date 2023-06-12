"""

This script labels CA to 1, N to 2 and C to 3

"""
import sys

import mrcfile
import math
import numpy as np
from Bio import PDB
import os

error_list = set()


def get_index(cord, origin, voxel):
    return math.ceil(math.floor(cord - origin) / voxel)

def ca_mask(path, org_map, pdb_map, outfilename):
    count = 0
    org_map = os.path.join(path, org_map)
    org_map = mrcfile.open(org_map, mode='r')
    data = np.zeros(org_map.data.shape, dtype=np.int16)
    data = data.astype('float32')
    x_origin = org_map.header.origin['x']
    y_origin = org_map.header.origin['y']
    z_origin = org_map.header.origin['z']
    x_voxel = org_map.voxel_size['x']
    y_voxel = org_map.voxel_size['y']
    z_voxel = org_map.voxel_size['z']
    parser = PDB.PDBParser()
    pdb_map = os.path.join(path, pdb_map)
    struct = parser.get_structure("CA", pdb_map)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x, y, z = atom.get_coord()
                    iz = int(get_index(z, z_origin, z_voxel))
                    jy = int(get_index(y, y_origin, y_voxel))
                    kx = int(get_index(x, x_origin, x_voxel))
                    if atom.get_name() == "CA":
                        try:
                            data[iz, jy, kx] = 1
                            count += 1
                        except IndexError as error:
                            error_list.add(pdb_map)
    print("Saving the file - - - - - ", org_map)
    outfilename = path + "/" + outfilename
    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = org_map.header.origin
        mrc.close()
    print(outfilename, "Done")
    print(f"Number of Carbon Alpha in label map is => {count}")

def label_mask(path, org_map, pdb_map, outfilename):
    count = 0
    org_map = os.path.join(path, org_map)
    org_map = mrcfile.open(org_map, mode='r')
    data = np.zeros(org_map.data.shape, dtype=np.int16)
    data = data.astype('float32')
    x_origin = org_map.header.origin['x']
    y_origin = org_map.header.origin['y']
    z_origin = org_map.header.origin['z']
    x_voxel = org_map.voxel_size['x']
    y_voxel = org_map.voxel_size['y']
    z_voxel = org_map.voxel_size['z']
    parser = PDB.PDBParser()
    pdb_map = os.path.join(path, pdb_map)
    struct = parser.get_structure("CA", pdb_map)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x, y, z = atom.get_coord()
                    iz = int(get_index(z, z_origin, z_voxel))
                    jy = int(get_index(y, y_origin, y_voxel))
                    kx = int(get_index(x, x_origin, x_voxel))
                    if atom.get_name() == "CA":
                        try:
                            data[iz, jy, kx] = 1
                            count += 1
                        except IndexError as error:
                            error_list.add(pdb_map)
                    elif atom.get_name() == "N":
                        try:
                            data[iz, jy, kx] = 2
                        except IndexError as error:
                            error_list.add(pdb_map)
                    elif atom.get_name() == "C":
                        try:
                            data[iz, jy, kx] = 3
                        except IndexError as error:
                            error_list.add(pdb_map)
    print("Saving the file - - - - - ", org_map)
    outfilename = path + "/" + outfilename
    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = org_map.header.origin
        mrc.close()
    print(outfilename, "Done")
    print(f"Number of Carbon Alpha in label map is => {count}")


if __name__ == "__main__":
    input_path = sys.argv[1]
    undone_pdb_emd = list()
    map_names = [fn for fn in os.listdir(input_path) if not fn.endswith(".pdb")]
    print("########### Generating atoms label ##########")
    for _ in range(len(map_names)):
        path = os.path.join(input_path, map_names[_])
        emd_map = [e for e in os.listdir(path) if e.endswith(".mrc")]
        pdb_map = [p for p in os.listdir(path) if p.endswith(".pdb")]
        pdb_map.sort()
        pdb_map = pdb_map[0].split(".")[0]
        pdb_map = pdb_map.split("_")[0]
        pdb_map = pdb_map.lower()
        pdb_map = pdb_map + ".pdb"
        if len(pdb_map) != 0 and len(emd_map) != 0:
            print("Working on => ", pdb_map, " of => ", map_names[_])
            em = "emd_normalized_map.mrc"
            name = em.split(".")
            # atoms labeling
            label_mask(path, em, pdb_map, "atom_" + name[0] + ".mrc")
            # ca only labeling
            ca_mask(path, em, pdb_map, "atom_ca_" + name[0] + ".mrc")
        else:
            undone_pdb_emd.append(map_names[_])

    print("Atoms and Ca only labeling Complete!")

