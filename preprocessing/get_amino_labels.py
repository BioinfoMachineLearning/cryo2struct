"""
@author: nabin

This script uses normalized emd map and corresponding pdb file to extract amino acid from pdb file and create a new
mrc file with amino acids.
For 20 difference amino acids, we label each one of them with numbers, 0 means no presence of amino acid, 
1 means ALA, 2 means ARG, 3 means ASP and so on. Which means we label total 21 -> 20 different amino acids and 0 means 
no presence of amino acids.

Note: GLY does not contain "CB" hence, we use "CA" for all.

"""
import sys
import mrcfile
import math
import numpy as np
from Bio import PDB
import os

error_set = set()

# A list of atoms (excluding hydrogen) for each AA type. PDB naming convention.
residue_atoms = {
    'ALA': ['C', 'CA', 'CB', 'N', 'O'],
    'ARG': ['C', 'CA', 'CB', 'CG', 'CD', 'CZ', 'N', 'NE', 'O', 'NH1', 'NH2'],
    'ASP': ['C', 'CA', 'CB', 'CG', 'N', 'O', 'OD1', 'OD2'],
    'ASN': ['C', 'CA', 'CB', 'CG', 'N', 'ND2', 'O', 'OD1'],
    'CYS': ['C', 'CA', 'CB', 'N', 'O', 'SG'],
    'GLU': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'O', 'OE1', 'OE2'],
    'GLN': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'NE2', 'O', 'OE1'],
    'GLY': ['C', 'CA', 'N', 'O'],
    'HIS': ['C', 'CA', 'CB', 'CG', 'CD2', 'CE1', 'N', 'ND1', 'NE2', 'O'],
    'ILE': ['C', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'O'],
    'LEU': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'N', 'O'],
    'LYS': ['C', 'CA', 'CB', 'CG', 'CD', 'CE', 'N', 'NZ', 'O'],
    'MET': ['C', 'CA', 'CB', 'CG', 'CE', 'N', 'O', 'SD'],
    'PHE': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'N', 'O'],
    'PRO': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'O'],
    'SER': ['C', 'CA', 'CB', 'N', 'O', 'OG'],
    'THR': ['C', 'CA', 'CB', 'CG2', 'N', 'O', 'OG1'],
    'TRP': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3',
            'CH2', 'N', 'NE1', 'O'],
    'TYR': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'N', 'O',
            'OH'],
    'VAL': ['C', 'CA', 'CB', 'CG1', 'CG2', 'N', 'O']
}
restype_1to3 = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}
residue_label = {
    'ALA': 1,
    'ARG': 2,
    'ASN': 3,
    'ASP': 4,
    'CYS': 5,
    'GLN': 6,
    'GLU': 7,
    'GLY': 8,
    'HIS': 9,
    'ILE': 10,
    'LEU': 11,
    'LYS': 12,
    'MET': 13,
    'PHE': 14,
    'PRO': 15,
    'SER': 16,
    'THR': 17,
    'TRP': 18,
    'TYR': 19,
    'VAL': 20,
}


def get_index(cord, origin, voxel):  # formula to convert coordinates to index
    return math.ceil(math.floor(cord - origin) / voxel)


def label_generator(path, org_map, pdb_map, outfilename):
    org_map = os.path.join(path, org_map)
    org_map = mrcfile.open(org_map, mode="r")
    data = np.zeros(org_map.data.shape, dtype=np.float32)
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
                    if atom.get_name() == "CA":
                        x, y, z = atom.get_coord()
                        iz = int(get_index(z, z_origin, z_voxel))
                        jy = int(get_index(y, y_origin, y_voxel))
                        kx = int(get_index(x, x_origin, x_voxel))
                        try:
                            label = residue_label[residue.resname]
                        except KeyError:
                            label = 0
                        try:
                            data[iz, jy, kx] = label
                        except IndexError as error:
                            error_set.add(pdb_map)
    outfilename = path + "/" + outfilename
    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = org_map.header.origin
        mrc.close()
    print(f"The voxel size is {x_voxel}")
    print(f"The saved filename is {outfilename}")


if __name__ == "__main__":
    count = 0
    undone_pdb_emd = list()
    input_path = sys.argv[1]  # path to normalized raw data - emd and pdb
    map_names = [protein for protein in os.listdir(input_path) if
                 not protein.endswith(".ent")]  # keep proteins in list except for the bounding box
    
    print("###### Extracting amino acid labels ######")
    for proteins in range(len(map_names)):
        path = os.path.join(input_path, map_names[proteins])
        try:
            emd = [e for e in os.listdir(path) if e.endswith(".mrc") and os.path.isdir(path)]
            pdb = [p for p in os.listdir(path) if p.endswith(".pdb") and os.path.isdir(path)]
        except NotADirectoryError:
            continue
        print(f"The corresponding PDB map for {emd} is {pdb}")
        print(path)

        if len(pdb) != 0 and len(emd) != 0:
            em = "emd_normalized_map.mrc"
            map_name = em.split(".")
            label_generator(path, em, pdb[0], "amino_" + map_name[0] + ".mrc")
            count += 1
        else:
            undone_pdb_emd.append(map_names[proteins])

    print("Label generation for amino acids complete!")
    print("Total done maps:", count)
    print('Undone Maps:', undone_pdb_emd)
