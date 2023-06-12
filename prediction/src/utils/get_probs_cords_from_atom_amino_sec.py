"""
Created on 02 March 2023 05:14:00 PM
@author: nabin

This script gets the carbon alphas from the amino prediction. Each voxel has probabilities of amino acids, so we sum the amino and if its higher than 0.5
then that means, that voxel has high chance of having carbon alpha atom and use that for HMM
"""
import ast
import os
import math

ca_coordinates = list()
prob_dic = dict()


def get_joint_probabity_common_sec_threshold(probability_file_atom, probability_file_amino_common_emi, probability_file_sec_common_emi, probability_file_amino, probability_file_sec, s_c, threshold):
    """
    get only common carbon alphas from atom, amino and sec files
    
    """
    common_ca = dict()
    common_coordinate_atm_sec = dict()
    amino_acid_emission = dict()
    sec_emission = dict()
    count_uncommon_atoms = 0
    total_sec_atom_common = 0
    total_saved_ca = 0 


    with open(probability_file_atom, 'r') as atom_prob:
        for line in atom_prob:
            line_a = ast.literal_eval(line)
            common_ca[tuple(line_a[0])] = line_a[2]


    with open(probability_file_sec, 'r') as sec_prob:
        for line in sec_prob:
            line_a = ast.literal_eval(line)
            ca_cord = tuple(line_a[0])
            ca_prob = 1 - line_a[1]
            try:
                test_cord =  math.sqrt(common_ca[ca_cord] * ca_prob)
                common_coordinate_atm_sec[ca_cord] = common_ca[ca_cord]
                ss_val = list(line_a[2:])
                sec_emission[ca_cord] = ss_val
                total_sec_atom_common += 1
            except KeyError:
                count_uncommon_atoms += 1
    
    with open(probability_file_amino, 'r') as amino_prob:
        for line in amino_prob:
            line_a = ast.literal_eval(line)
            ca_cord = tuple(line_a[0])
            ca_prob = 1 - line_a[1]
            try:
                test_cord = math.sqrt(common_coordinate_atm_sec[ca_cord] * ca_prob)
                aa_val = list(line_a[2:])
                amino_acid_emission[ca_cord] = aa_val
                # total_saved_ca += 1
            except KeyError:
                count_uncommon_atoms += 1


    print("Using threshold for predicted probabilities ==>", threshold)
    save_cluster_co = open(s_c, 'a')
    amino_emi = open(probability_file_amino_common_emi,'a')
    sec_emi = open(probability_file_sec_common_emi, 'a')
    for k,v in common_coordinate_atm_sec.items():
        if v > threshold:
            try:
                emiss_val = amino_acid_emission[k]
                emiss_val_ss = sec_emission[k]
                x,y,z = k
                print(f"{x} {y} {z}", file=save_cluster_co)
                amino_emi.write(f"{list(k)}")
                for e in emiss_val:
                    amino_emi.write(f", {e}")
                amino_emi.write(f"\n")

                sec_emi.write(f"{list(k)}")
                for s in emiss_val_ss:
                    sec_emi.write(f", {s}")
                sec_emi.write(f"\n")

                total_saved_ca += 1
            except KeyError:
                q  = 1

    print(f" Number of SAVED carbon alpha in {s_c} file ATOM AMINO SEC ==> {total_saved_ca}")
    print(f"Coordinate file saved in ==> ", s_c)
    print(f"Amino Emission file saved in ==> ", probability_file_amino_common_emi)
    print(f"Sec Emission file saved in ==> ", probability_file_sec_common_emi)

def get_joint_probabity_common_threshold(probability_file_atom, probability_file_amino_atom_common, probability_file_amino, s_c, threshold):
    """
    get only common carbon alphas from amino and atom files
    
    """
    common_ca = dict()
    common_coordinate_prob = dict()
    amino_acid_emission = dict()
    count_uncommon_atoms = 0
    count_common_atoms = 0
    total_atom_entries = 0
    total_saved_ca = 0

    with open(probability_file_amino, 'r') as amino_prob:
        for line in amino_prob:
            line_a = ast.literal_eval(line)
            common_coordinate_prob[tuple(line_a[0])] = 1 - line_a[1]
            aa_val = list(line_a[2:])
            equal_part_add = line_a[1]  / 20
            aa_val = tuple([x + equal_part_add for x in aa_val])
            amino_acid_emission[tuple(line_a[0])] = aa_val
            total_atom_entries += 1
                
            
    with open(probability_file_atom, 'r') as atom_prob:
        for line in atom_prob:
            line_a = ast.literal_eval(line)
            try:
                common_ca[tuple(line_a[0])] = math.sqrt(common_coordinate_prob[tuple(line_a[0])] * line_a[2])
                common_ca[tuple(line_a[0])] = line_a[2]
                count_common_atoms += 1 
    
            except KeyError:
                count_uncommon_atoms += 1

    print("Using threshold for predicted probabilities ==>", threshold)
    save_cluster_co = open(s_c, 'a')
    amino_atom_prob = open(probability_file_amino_atom_common,'a')
    for k,v in common_ca.items():
        if v > threshold:
            try:
                emiss_val = amino_acid_emission[k]
                x,y,z = k
                print(f"{x} {y} {z}", file=save_cluster_co)
                amino_atom_prob.write(f"{list(k)}")
                for e in emiss_val:
                    amino_atom_prob.write(f", {e}")
                amino_atom_prob.write(f"\n")
                total_saved_ca += 1
            except KeyError:
                q  = 1

    print(f" Number of TOTAL carbon alpha in amino file ==> {total_atom_entries}")
    print(f" Number of COMMON carbon alpha in atom and amino file ==> {count_common_atoms}")
    print(f" Number of SAVED carbon alpha in {s_c} file ==> {total_saved_ca}")
    print(f" Number of UNCOMMON carbon alpha in atom and amino file ==> {count_uncommon_atoms}")

