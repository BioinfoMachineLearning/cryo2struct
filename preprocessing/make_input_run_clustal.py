"""
@author: nabin

1 - This script reads combined fasta and atom fasta then, writes into dealign_clustal_input.fasta
2 - Runs clustal omega on dealign_clustal_input.fasta. Output -> dealign_clustal_output.fasta

>original
# original fasta seq here, <map_name>_all_chain_combined.fasta
>atom
# atom fasta seq here, atomic.fasta

"""
import os
import sys


# store the density map names having issues:
string_to_check_found_density_maps = set()


def make_clustal_input(em_maps, input_path, string_to_check):
    print("Making clustal inputs!!!")
    chain_fastas = [fa for fa in os.listdir(os.path.join(input_path, em_maps)) if fa.endswith(".fasta")]
    chain_fastas.sort()
    f_name = chain_fastas[0].split(".")[0]
    header_original = ">original"
    header_atom = ">atom"
    original_fasta = f"{f_name}_all_chain_combined.fasta"
    atom_fasta = "atomic.fasta"
    filename = f'{input_path}/{em_maps}/dealign_clustal_input.fasta'
    file_open = open(filename, "a")
    file_open.write(header_original)
    file_open.write("\n")
    with open(os.path.join(input_path, em_maps, original_fasta)) as org:
        lines = org.readlines()
    file_open.write(lines[0])
    file_open.write("\n")
    file_open.write(header_atom)
    file_open.write("\n")
    with open(os.path.join(input_path, em_maps, atom_fasta)) as atom:
        lines = atom.readlines()
    with open(os.path.join(input_path, em_maps, atom_fasta), 'w') as file:
        for l in lines:
            # Check if the line starts with the specified string
            if not l.startswith(string_to_check):
                # If not, write the line to the file
                file.write(l)
            else:
                string_to_check_found_density_maps.add(os.path.join(input_path, em_maps))
    with open(os.path.join(input_path, em_maps, atom_fasta)) as atom:
        lines = atom.readlines()
    file_open.write(lines[0])
    file_open.close()
    print(filename, "-> Done")


def run_clustal(em_maps, input_path):
    print("Running Clustal Now!!! ")
    density_maps_dir = os.path.join(input_path, em_maps)
    input_clustal_file = os.path.join(density_maps_dir, "dealign_clustal_input.fasta")
    os.system(
        "clustalo -i " + input_clustal_file + " --dealign" + ">>" + density_maps_dir + "/" + "dealign_clustal_output.fasta")
    print("Done ->", em_maps)


if __name__ == "__main__":
    input_path = sys.argv[1]
    # remove if any of the atom fasta starts with the below:
    string_to_check = f"{input_path}/"
    count = 0
    density_maps = [d_map for d_map in os.listdir(input_path) if
                    os.path.isdir(os.path.join(input_path, d_map))]

    for em_maps in density_maps:
        rm1 = f'{input_path}/{em_maps}/dealign_clustal_input.fasta'
        if os.path.exists(rm1):
            os.remove(rm1)
        rm2 = f'{input_path}/{em_maps}/dealign_clustal_output.fasta'
        if os.path.exists(rm2):
            os.remove(rm2)
        make_clustal_input(em_maps, input_path, string_to_check)
        run_clustal(em_maps, input_path)
        print(count)
        count += 1
    print(string_to_check_found_density_maps)
