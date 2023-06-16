"""
@author: nabin

This script extracts helix, strand and coil from the pdb file and saves them separately

"""
import subprocess
import os
import sys

secondary_atoms = ["helix", "strand", "coil"]

def extract(input_maps, maps, chimera_path):
    density_map = os.path.join(input_maps, maps)
    pdbs = [pdb for pdb in os.listdir(density_map) if pdb.endswith(".pdb")]
    print(f"The pdb files found are : {pdbs}")
    if pdbs is None:
        print("Please check the directory and rerun, there are no pdb files found")
        exit()
    else:
        for pdb in range(len(pdbs)):
            for sec_atom in secondary_atoms:
                chimera_scripts = open('sec_extractor.cxc', 'w')
                chimera_scripts.write('open ' + density_map + '/' + pdbs[pdb] + '\n'
                                      'select ' + sec_atom + '\n'
                                      'save ' + density_map + "/" + sec_atom +'.pdb' + ' format pdb ' + 'selectedOnly true \n'
                                      'exit')
                chimera_scripts.close()
                script_finished = False
                while not script_finished:
                    try:
                        subprocess.run([chimera_path,'--nogui', chimera_scripts.name])
                        script_finished = True
                    except FileNotFoundError as error:
                        raise error
                os.remove(chimera_scripts.name)
                print(f"Done {sec_atom} for {pdbs[pdb]}")
        print("NG: All secondary atom PDBs are Done")



if __name__ == "__main__":
    input_maps = sys.argv[1]
    if len(sys.argv) > 2:
        chimera_path = sys.argv[2]
    else:
        chimera_path = '/usr/bin/chimerax'

    density_maps = os.listdir(input_maps)
    print(f"Found density maps are :{density_maps}")
    for maps in density_maps:
        print(f"Working on {maps}")
        extract(input_maps, maps, chimera_path)
