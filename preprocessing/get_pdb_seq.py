"""
@author: nabin

This script runs pdb2seq.pl to get sequence from pdb file

"""
import os
import sys


def run_pdb2sec(pdb, density_name, input_path, perl_script_dir):
    try:
        density_map_dir = os.path.join(input_path, density_name)
        pdb_fi = os.path.join(density_map_dir, pdb)
        os.system("perl " + perl_script_dir + " " + pdb_fi + ">>" + density_map_dir + "/" + "atomic.fasta")
        print(density_name, "Done")
    except FileNotFoundError:
        print("Error for file:", density_name)


if __name__ == "__main__":
    input_path = sys.argv[1]
    perl_script_dir = "pdb2seq.pl"
    density_maps = [den for den in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, den))]
    for den in density_maps:
        rm1 = f'{input_path}/{den}/atomic.fasta'
        if os.path.exists(rm1):
            os.remove(rm1)
        pdb_file = [p for p in os.listdir(os.path.join(input_path, den)) if
                    p.endswith(".pdb") or p.endswith(".ent")]
        pdb_file.sort()
        pdb_file_n = pdb_file[0].split(".")[0]
        pdb_file_na = pdb_file_n.split("_")[0]
        pdb_file_name = pdb_file_na + ".pdb"
        # print(den, "->", pdb_file_name)

        run_pdb2sec(pdb_file_name, den, input_path, perl_script_dir)
