"""
Created on 18 April 2023 12:23 AM
@author: nabin

Usage:
- Gets sequence from pdb file along with its chain information

"""


from Bio import PDB

restype_3to1 = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'UNK' : 'U',
}


chain_seq_dict = dict()

def extract_seq(pdb_file, atomic_chain_seq_file, atomic_seq_file):
    parser = PDB.PDBParser()
    pdb_map = pdb_file
    struct = parser.get_structure("CA", pdb_map)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "CA":
                        chain_id = chain.id
                        try:
                            amino_name = restype_3to1[residue.resname]
                            if chain_id in chain_seq_dict:
                                chain_seq_dict[chain_id].append(amino_name)
                            else:
                                chain_seq_dict[chain_id] = [amino_name]
                        except KeyError:
                            pass
                        
    with open(atomic_chain_seq_file, 'w') as a_c:
        for k,v in chain_seq_dict.items():
            print(f">pdb2seq|Chains {k}", file=a_c)
            result = ''.join(v)
            print(result, file=a_c)
    
    all_seq = list()
    with open(atomic_seq_file, 'w') as a_s:
        print(">pdb2seq|Chains A", file=a_s)
        for k,v in chain_seq_dict.items():
            result = ''.join(v)
            all_seq.append(result)
        final_result = ''.join(all_seq)
        print(final_result,file=a_s)
