"""
Gets FASTA file of the proteins from RCSB Website.
"""

import urllib.request
import urllib.error
import os
import requests

def get_pdb(em_path):
    counter = 0
    dir_names = [m for m in os.listdir(em_path)]
    print("Length of maps", len(dir_names))
    for e in dir_names:
        density = os.path.join(em_path,e)
        pdb_files = [p for p in os.listdir(density) if p.endswith(".pdb")]
        pdb_files.sort()
        pdb_files = pdb_files[0].split("_")
        pdb_id = pdb_files[0].split(".")[0]
        
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"

        # Send an HTTP GET request to the URL
        response = requests.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            # Extract the content of the response
            fasta_content = response.text

            # Save the content to a file
            filename = f"{density}/{pdb_id}.fasta"
            with open(filename, "w") as file:
                file.write(fasta_content)
            counter += 1

            # print(f"FASTA file saved as: {filename}")
        else:
            print(f"Failed to download FASTA file for PDB ID: {pdb_id} => {density}")


    print("ALL DONE : ", counter)




if __name__ == "__main__":
    emd_path = sys.argv[1]
    get_pdb(emd_path)