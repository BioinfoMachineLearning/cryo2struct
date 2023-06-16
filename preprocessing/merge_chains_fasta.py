"""
@author: nabin

- This script merges the chains of original fasta to a single sequence and saves it.
- Use chain_merger if you have chain sequence on multiple fasta files (<pdb_code>_1.fasta, <pdb_code>_2.fasta, <pdb_code>_3.fasta and so on).
- Use chain_merger_2 if you have chain sequences on a single fasta file (<pdb_code>.fasta).

"""
import os
import sys


def chain_merger_2(map_name, input_path):
    chain_fastas = [fa for fa in os.listdir(os.path.join(input_path, map_name)) if fa.endswith(".fasta")]
    chain_fastas.sort()
    f_name = chain_fastas[0].split("_")
    f_name = chain_fastas[0].split(".")[0]
    input_file = os.path.join(input_path,map_name,f'{f_name}.fasta')


    output_file = f'{input_path}/{map_name}/{f_name}_all_chain_combined.fasta'
    if os.path.exists(output_file):
        os.remove(output_file)
    
    with open(input_file, "r") as input_fp, open(output_file, "w") as output_fp:
        merge_lines = []
        for line in input_fp:
            if line.startswith(">"):
                if merge_lines:
                    merged_line = "".join(merge_lines)
                    output_fp.write(merged_line)
                    merge_lines = []
            else:
                merge_lines.append(line.strip())
        if merge_lines:
            merged_line = "".join(merge_lines)
            output_fp.write(merged_line)
    print(f_name, "Done")

def save_combined_fasta(map_name, combined_seq, input_path,f_name):
    filename = f'{input_path}/{map_name}/{f_name}_all_chain_combined.fasta'
    file_open = open(filename, "w")
    file_open.write(combined_seq)
    file_open.close()
    print(filename, "-> Done")



def chain_merger(map_name, input_path):
    """
    :type map_name: density map name

    create a dictionary with chain name as key and chain sequence as values
    """
    chain_seq = dict()
    chain_fastas = [fa for fa in os.listdir(os.path.join(input_path, map_name)) if fa.endswith(".fasta")]
    chain_fastas.sort()
    f_name = chain_fastas[0].split(".")[0]
    for fa in chain_fastas:
        print(fa)
        with open(os.path.join(input_path, map_name, fa)) as f:
            lines = f.readlines()
            header = lines[0]
            seq = lines[1]
            chains_ids = header.split("|")[1]
            chains = chains_ids.split(",")
            for ch in chains:
                ch = ch.replace("Chains", "")
                ch = ch.replace("Chain", "")
                ch = ch.strip(' ')
                seq = seq.rstrip('\n')
                chain_seq[ch] = seq
    chain_seq = sorted(chain_seq.items())
    chain_seq = dict(chain_seq)

    combined_seq = ""
    length = 0
    for value in chain_seq.values():
        combined_seq += value
        length += len(value)
    save_combined_fasta(map_name, combined_seq, input_path, f_name)
    print("Validate the length of sequences: ", length == len(combined_seq))


if __name__ == "__main__":
    input_path = sys.argv[1]
    density_maps = [d_map for d_map in os.listdir(input_path) if
                    os.path.isdir(os.path.join(input_path, d_map))]
    for em_maps in density_maps:
        chain_merger_2(em_maps, input_path)
    print("Chain merger completed!")
