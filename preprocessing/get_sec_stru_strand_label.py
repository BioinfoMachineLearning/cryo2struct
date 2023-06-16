"""
@author: nabin

This script opens in read/write mode of sec_struc_emd_normalized_map.mrc MRC file and labels strand to 3
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


def label_mask(path, org_map, pdb_map, outfilename):
    org_map = os.path.join(path, org_map)
    org_map = mrcfile.open(org_map, mode='r')
    x_origin = org_map.header.origin['x']
    y_origin = org_map.header.origin['y']
    z_origin = org_map.header.origin['z']
    x_voxel = org_map.voxel_size['x']
    y_voxel = org_map.voxel_size['y']
    z_voxel = org_map.voxel_size['z']
    parser = PDB.PDBParser()
    pdb_map = os.path.join(path, pdb_map)
    struct = parser.get_structure("CA", pdb_map)
    outfilename = path + "/" + outfilename
    with mrcfile.open(outfilename, mode='r+') as mrc:
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
                                mrc.data[iz, jy, kx] = 3
                            except IndexError as error:
                                error_list.add(pdb_map)

        print("Saving the file - - - - - ", org_map)
        mrc.close()
    print(outfilename, "Done")


if __name__ == "__main__":
    input_path = sys.argv[1]
    undone_pdb_emd = list()
    count = 0
    map_names = [fn for fn in os.listdir(input_path) if not fn.endswith(".DS_Store")]
    print("########### Generating SS-strand labels ##########")
    for _ in range(len(map_names)):
        path = os.path.join(input_path, map_names[_])
        emd_map = [e for e in os.listdir(path) if e.endswith(".mrc")]
        pdb_map = [p for p in os.listdir(path) if p.endswith(".pdb")]
        if len(pdb_map) != 0 and len(emd_map) != 0:
            em = "emd_normalized_map.mrc"
            name = em.split(".")
            pdb_map = pdb_map[0].split(".")
            pdb_map = pdb_map[0].split("_")
            pdb_map = pdb_map[0]
            pdb_map = f"strand.pdb"
            label_mask(path, em, pdb_map, "sec_struc_" + name[0] + ".mrc")
            count += 1
        else:
            undone_pdb_emd.append(map_names[_])

    print("Secondary structure labeling for : Strands - Complete!")
    print("Total done maps:", count)
    print("Undone maps:", undone_pdb_emd)
