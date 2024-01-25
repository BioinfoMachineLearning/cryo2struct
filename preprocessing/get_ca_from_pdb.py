"""

@author: nabin

This script extracts carbon-alpha only atoms from protein structure (.pdb) and saves it.
Run:
python3 get_ca_from_pdb.py <input protein structure path> <output protein only ca structure path>

"""

from Bio import PDB
import os
import sys

def save(x, y, z, count, amino_name, chain_id, true_ca_path):
    atom = 'CA'
    residue_name = amino_name
    with open(true_ca_path, 'a') as fi:
        fi.write('ATOM')
        fi.write('  ')
        fi.write(str(count).rjust(5))
        fi.write('  ')
        fi.write(atom.ljust(4))
        fi.write(residue_name.rjust(3))
        fi.write(' ')
        fi.write(chain_id)  
        fi.write(str(count).rjust(4))
        fi.write('    ')
        fi.write(str(x).rjust(8))
        fi.write(str(y).rjust(8))
        fi.write(str(z).rjust(8))
        fi.write(str(1.00).rjust(5))
        fi.write(str(0.00).rjust(5))
        fi.write('           ')
        fi.write(atom[0:1].rjust(1))
        fi.write('  ')
        fi.write('\n')


def label_mask(pdb_path, true_ca_path):
    count = 0
    parser = PDB.PDBParser(QUIET=True)
    pdb_map = pdb_path
    struct = parser.get_structure("CA", pdb_map)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "CA":
                        x, y, z = atom.get_coord()
                        amino_name = residue.resname
                        chain_id = chain.id
                        save(x, y, z, count, amino_name, chain_id, true_ca_path)
                        count += 1


if __name__ == "__main__":
    input_pdb = sys.argv[1]
    true_ca_pdb = sys.argv[2]
    if os.path.exists(true_ca_pdb):
            os.remove(true_ca_pdb)
    label_mask(input_pdb, true_ca_pdb)
    print("########## NG: get_ca_labels: DONE ##########")

