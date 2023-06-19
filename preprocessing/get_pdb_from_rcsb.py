
"""
@author: nabin

Gets PDB file of the proteins from RCSB Website.
"""

import urllib.request
import urllib.error
import os
import sys
import requests


def get_pdb(em_path):
    counter = 0
    dir_names = [m for m in os.listdir(em_path)]
    print("Length of maps", len(dir_names))
    for e in dir_names:
        # Retrieve the PDB structures associated with the EMDB entry
        api_url = f"https://www.ebi.ac.uk/emdb/api/entry/{e}"
        response = requests.get(api_url)

        if response.status_code == 200:
            data = response.json()
            if "crossreferences" in data and "pdb_list" in data["crossreferences"]:
                pdb_entries = data["crossreferences"]["pdb_list"]
            
                pdb_id = pdb_entries['pdb_reference']
                pdb_id = pdb_id[0]['pdb_id']

                # Download the PDB file
                pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                pdb_response = requests.get(pdb_url)

                if pdb_response.status_code == 200:
                    pdb_filename = f"{os.path.join(em_path,e)}/{pdb_id}.pdb"
                    with open(pdb_filename, "wb") as file:
                        file.write(pdb_response.content)
                    print(f"Successfully downloaded PDB file: {pdb_filename}")
                else:
                    print(f"Error downloading the PDB file for PDB ID: {pdb_id}")
        else:
            print("No corresponding PDB structures found for the given EMDB ID.")
    else:
        print("Error retrieving data from the EMDB API.")


    print("DONE : ", counter)




if __name__ == "__main__":
    emd_path = sys.argv[1]
    get_pdb(emd_path)
