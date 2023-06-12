"""
Created on 17 April 2023 10:47 AM
@author: nabin

Usage:
- Connects carbon-alpha atom to from atomic backbone structure for protein

"""


import argparse
import yaml
import os
import mrcfile


from infer import ca_final_amino_inference_no_attn, ca_final_atom_inference_no_attn 
from utils import get_ca_from_pred_probs, get_probs_cords_from_atom_amino_sec, clustering_centroid
from viterbi import viterbi_batches


script_dir = os.path.dirname(os.path.abspath(__file__))

config_file_path = f"{script_dir}/config/arguments.yml"
COMMENT_MARKER = '#'

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=argparse.FileType(mode='r'),
                        default=config_file_path)
    parser.add_argument('--density_map_name', type=str)
    
    return parser.parse_args()


def process_arguments(args):
    if args.config is not None:
        config_dict = yaml.safe_load(args.config)
        config_dict = {k: v for k, v in config_dict.items() if not k.startswith(COMMENT_MARKER)}
        args.config = args.config.name
    else:
        config_dict = dict()
    
    if args.density_map_name is not None:
        config_dict['density_map_name'] = args.density_map_name

    print(config_dict)

    return config_dict
    
    
def make_predictions(config_dict):
    print("##############- MAKING PREDICTIONS -##############")
    ca_final_amino_inference_no_attn.infer_node_classifier(test_data_splits_dir=config_dict['test_data_splits_dir'], 
                                                           test_data_dir=config_dict['test_data_dir'], 
                                                           density_map_name=config_dict['density_map_name'], 
                                                           amino_checkpoint=config_dict['amino_checkpoint'], 
                                                           infer_run_on=config_dict['infer_run_on'],
                                                           infer_on_gpu=config_dict['infer_on_gpu'])
    
    ca_final_atom_inference_no_attn.infer_node_classifier(test_data_splits_dir=config_dict['test_data_splits_dir'], 
                                                        test_data_dir=config_dict['test_data_dir'], 
                                                        density_map_name=config_dict['density_map_name'], 
                                                        atom_checkpoint=config_dict['atom_checkpoint'], 
                                                        infer_run_on=config_dict['infer_run_on'],
                                                        infer_on_gpu=config_dict['infer_on_gpu'])
    print("##############- MAKING PREDICTIONS COMPLETED -##############")

def extract_ca_from_prediction_probabilities(config_dict):
    original_map = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/emd_normalized_map.mrc"
    original_map_mrc = mrcfile.open(original_map, mode='r')
    original_map_shape = original_map_mrc.data.shape
    original_map_origin = original_map_mrc.header.origin
    x_origin = original_map_mrc.header.origin['x']
    y_origin = original_map_mrc.header.origin['y']
    z_origin = original_map_mrc.header.origin['z']
    x_voxel = original_map_mrc.voxel_size['x']
    y_voxel = original_map_mrc.voxel_size['y']
    z_voxel = original_map_mrc.voxel_size['z']


    # extract ca from atoms prediction : for visualization and evaluation
    pred_atom_prob = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_atom.txt"
    only_ca_atom_mrc = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_atom_ca_predicted.mrc"
    only_ca_atom_pdb = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_atom_ca_predicted.pdb"
    if os.path.exists(only_ca_atom_pdb):
        os.remove(only_ca_atom_pdb)
    get_ca_from_pred_probs.extract_ca_from_atom(pred_atom=pred_atom_prob, outfilename=only_ca_atom_mrc, outfilename_pdb=only_ca_atom_pdb, density_shape=original_map_shape, density_voxel=(x_voxel,y_voxel,z_voxel), density_origin=(x_origin,y_origin,z_origin), origin=original_map_origin)
    print("Extracting carbon alpha from ATOMS prediction complete!")


    # extract ca from amino prediction: for visulalization, and evaluation
    pred_amino_prob = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_amino.txt"
    only_ca_amino_mrc = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_amino_ca_predicted.mrc"
    only_ca_amino_pdb = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_amino_ca_predicted.pdb"
    if os.path.exists(only_ca_amino_pdb):
        os.remove(only_ca_amino_pdb)
    get_ca_from_pred_probs.extract_ca_from_amino(pred_amino=pred_amino_prob, outfilename=only_ca_amino_mrc, outfilename_pdb=only_ca_amino_pdb, density_shape=original_map_shape, density_voxel=(x_voxel,y_voxel,z_voxel), density_origin=(x_origin,y_origin,z_origin), origin=original_map_origin)
    print("Extracting carbon alpha from AMINO prediction complete!")


    # extract common ca from atoms and amino prediction : for visualization, and evaluation
    only_ca_common_from_atom_amino_mrc = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_atom_amino_ca_predicted.mrc"
    only_ca_common_from_atom_amino_pdb = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_atom_amino_ca_predicted.pdb"
    if os.path.exists(only_ca_common_from_atom_amino_pdb):
        os.remove(only_ca_common_from_atom_amino_pdb)
    get_ca_from_pred_probs.extract_ca_from_atom_amino_common_only(pred_atom=pred_atom_prob, pred_amino=pred_amino_prob, outfilename=only_ca_common_from_atom_amino_mrc, 
                                           outfilename_pdb=only_ca_common_from_atom_amino_pdb,
                                           density_shape=original_map_shape, density_voxel=(x_voxel,y_voxel,z_voxel), 
                                           density_origin=(x_origin,y_origin,z_origin), origin=original_map_origin)
    print("Extracting carbon alpha from ATOM AND AMINO prediction complete!")


    # extract common ca from atoms, amino and sec prediction : for visualization, and evaluation
    pred_sec_prob = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_sec.txt"
    only_ca_common_from_atom_amino_sec_mrc = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_atom_amino_sec_ca_predicted.mrc"
    only_ca_common_from_atom_amino_sec_pdb = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_prob_atom_amino_sec_ca_predicted.pdb"
    if os.path.exists(only_ca_common_from_atom_amino_sec_pdb):
        os.remove(only_ca_common_from_atom_amino_sec_pdb)
    print("Extracting carbon alpha from ATOM, AMINO AND SEC prediction complete!")


def extract_probs_cords_from_atom_amino_sec(config_dict):
    probability_file_atom = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_atom.txt" # comes from ca_final_atom_inference.py
    probability_file_amino = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_amino.txt" # comes from ca_final_amino_inference.py
    probability_file_sec = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_sec.txt" # comes from ca_final_sec_inference.py
    probability_file_amino_atom_common_emi = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_amino_atom_common_emi.txt" # save common amino and atom
    probability_file_amino_atom_sec_common_emi = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_amino_atom_sec_common_emi.txt" # save common amino, atom, and sec
    probability_file_amino_common_emi = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_amino_emi.txt" # save amino probability as emission
    probability_file_sec_common_emi = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_sec_emi.txt" # save sec probability as emission
    
    save_cords = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_coordinates_ca.txt" # save cords as transition matrix


    if os.path.exists(save_cords):
        os.remove(save_cords)
    
    if os.path.exists(probability_file_amino_atom_common_emi):
        os.remove(probability_file_amino_atom_common_emi)

    if os.path.exists(probability_file_amino_common_emi):
        os.remove(probability_file_amino_common_emi)

    if os.path.exists(probability_file_sec_common_emi):
        os.remove(probability_file_sec_common_emi)
        
    
    get_probs_cords_from_atom_amino_sec.get_joint_probabity_common_threshold(probability_file_atom=probability_file_atom, probability_file_amino_atom_common=probability_file_amino_atom_common_emi, 
                    probability_file_amino=probability_file_amino, s_c=save_cords, threshold = config_dict['threshold'])
    print("Saved in: ", probability_file_amino_atom_common_emi)


def cluster_emission_transition(config_dict):
    save_cords, save_probs_aa = clustering_centroid.main(config_dict)
    return save_cords, save_probs_aa


def main():
    args = parse_arguments()
    config_dict = process_arguments(args)

    if os.path.exists(f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}"):

        print("##############- RUNNING cryo2struct -##############")
        make_predictions(config_dict)

        # preparing for visualization and evaluation
        extract_ca_from_prediction_probabilities(config_dict)

        # preparing for HMM model 
        extract_probs_cords_from_atom_amino_sec(config_dict)

        # clustering and preparing emission and transition matrix
        coordinate_file, emission_file = cluster_emission_transition(config_dict)

        # run viterbi algorithm
        viterbi_batches.main(coordinate_file, emission_file, config_dict)
        print("ALL DONE")
    else:
        print("Missing density: ", config_dict['density_map_name'])
        exit()


if __name__ == "__main__":
    main()
