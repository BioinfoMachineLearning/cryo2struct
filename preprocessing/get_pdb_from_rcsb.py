
"""
@author: nabin

Gets PDB file of the proteins from RCSB Website.
"""

import urllib.request
import urllib.error
import os
import sys

def get_pdb(em_path):
    counter = 0
    dir_names = ['23552', '7635', '22883', '25263', '27068', '31760']
    print("Length of maps", len(dir_names))
    for e in dir_names:
        density = os.path.join(em_path,e)
        pdb_files = [p for p in os.listdir(density) if p.endswith(".pdb")]
        pdb_files.sort()
        pdb_files = pdb_files[0].split("_")
        pdb_files = pdb_files[0].split(".")[0]
        try:
            # Download the PDB file for the corresponding structure (PDB ID 6R4P)
            pdb_url = f'https://files.rcsb.org/download/{pdb_files}.fasta'
            urllib.request.urlretrieve(pdb_url, f'{density}/{pdb_files.lower()}.pdb')

            print(f'Successfully downloaded PDB file for EMD-ID {e} ===> {pdb_files.lower()}')
            counter += 1
            
        except urllib.error.HTTPError as e:
            print(f'Error downloading PDB file for EMD-ID {e} ===> {pdb_files}')

    print("DONE : ", counter)




if __name__ == "__main__":
    emd_path = sys.argv[1]
    get_pdb(emd_path)
