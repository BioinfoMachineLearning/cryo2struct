"""
Created on 18 April 2023 1:23 AM
@author: nabin

"""
import math
import ast
import mrcfile
import numpy as np
import os


prob_dic_aa = dict()
prob_dic_sec = dict()


class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def distance(p1, p2):
    dis =  math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)
    return dis


def create_clusters(points, thres):
    clusters = []
    while points:
        # Select a point randomly and assign it to a new cluster
        cluster = [points.pop(0)]
        # Iterate through the rest of the points and add them to the cluster if they are within the threshold distance
        i = 0
        while i < len(points):
            if distance(cluster[0], points[i]) <= thres:
                cluster.append(points.pop(i))
            else:
                i += 1
        clusters.append(cluster)
    return clusters


def centroid(file, save_cords, save_probs_aa, save_probs_sec , thres):
    # Read the data from the file
    points = []
    with open(file, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()

    for line in lines:
        vals = line.split(" ")
        for limiter in vals:
            if limiter == '':
                vals.remove(limiter)
        vals = list(filter(lambda x: x != '', vals))               
        points.append(Point(float(vals[0]), float(vals[1]), float(vals[2])))

    # Create the clusters
    clusters = create_clusters(points, thres=thres)
    print("Total coordinates => ", len(lines))
    print("Total clusters created => ", len(clusters))

    with open(save_probs_sec,'w') as s:
        with open(save_probs_aa,'w') as p:
            with open(save_cords, 'w') as f:
                for i, cluster in enumerate(clusters):
                    x_sum = 0
                    y_sum = 0
                    z_sum = 0
                    num_points = len(cluster)
                    collect_values = list()
                    collect_values_sec = list()
                    for point in cluster:
                        x_sum += point.x
                        y_sum += point.y
                        z_sum += point.z
                        cords = (point.x, point.y, point.z)
                        if cords in prob_dic_aa:
                            values = prob_dic_aa.get(cords)
                            collect_values.append(values)
                        if cords in prob_dic_sec:
                            values = prob_dic_sec.get(cords)
                            collect_values_sec.append(values)
                    averages = list()
                    averages_sec = list()
                    for i in range(len(collect_values[0])):
                        total = sum(collect_values[j][i] for j in range(len(collect_values)))
                        average = total / len(collect_values)     
                        averages.append(average) 
                    averages = ' '.join(str(x) for x in averages)
                    print(averages, file=p)   
                    x_avg = x_sum / num_points
                    y_avg = y_sum / num_points
                    z_avg = z_sum / num_points
                    print(f'{x_avg} {y_avg} {z_avg}', file=f)
                    
    
def proc_probabilities_sec(file):
    with open(file, 'r') as f:
        # Read all the lines in the file
        line = f.readline()
        while line:
            line_c = ast.literal_eval(line)
            key = tuple(line_c[0])
            vals = line_c[1:]
            prob_dic_sec[key] = vals
            line = f.readline()


def proc_probabilities_aa(file):
    with open(file, 'r') as f:
        # Read all the lines in the file
        line = f.readline()
        while line:
            line_c = ast.literal_eval(line)
            key = tuple(line_c[0])
            vals = line_c[1:]
            prob_dic_aa[key] = vals
            line = f.readline()

def save(x, y, z, count, save_path):
    atom = 'CA'
    residue_name = "GLY"
    with open(save_path, 'a') as fi:
        fi.write('ATOM')
        fi.write('  ')
        fi.write(str(count).rjust(5))
        fi.write('  ')
        fi.write(atom.ljust(4))
        fi.write(residue_name.rjust(3))
        fi.write(' ')
        fi.write('A') 
        fi.write(str(count).rjust(4))
        fi.write('    ')
        fi.write(str(x).rjust(8))
        fi.write(str(y).rjust(8))
        fi.write(str(z).rjust(8))
        fi.write(str(1.00).rjust(5))
        fi.write(str(0.00).rjust(5))
        fi.write('           ')
        fi.write(atom[0:1].rjust(1))
        fi.write('  ')
        fi.write('\n')


def get_index(cord, origin, voxel):
    return math.ceil(math.floor(cord - abs(origin)) / voxel)

def generate_pdb_mrc(save_cords, save_cords_pdb, save_cords_mrc, org_mrc):
    count = 0
    org_emd_mrc = mrcfile.open(org_mrc, mode="r")
    org_emd_mrc_data_shape = org_emd_mrc.data.shape
    origin_header = org_emd_mrc.header.origin
    x_origin = org_emd_mrc.header.origin['x']
    y_origin = org_emd_mrc.header.origin['y']
    z_origin = org_emd_mrc.header.origin['z']
    x_voxel = org_emd_mrc.voxel_size['x']
    y_voxel = org_emd_mrc.voxel_size['y']
    z_voxel = org_emd_mrc.voxel_size['z']
    data = np.zeros(org_emd_mrc_data_shape, dtype=np.float32)
    with open(save_cords,"r") as cluster_cords:
        line = cluster_cords.readlines()
    for l in line:
        count += 1
        vals = l.strip()
        vals = vals.split(" ")
        vals_x = round(float(vals[0]),3)
        vals_y = round(float(vals[1]),3)
        vals_z = round(float(vals[2]),3)
        save(x=vals_x, y=vals_y, z=vals_z, count=count, save_path=save_cords_pdb)
        iz = int(get_index(vals_z, z_origin, z_voxel))
        jy = int(get_index(vals_y, y_origin, y_voxel))
        kx = int(get_index(vals_x, x_origin, x_voxel))
        data[iz,jy,kx] = 1


    with mrcfile.new(save_cords_mrc, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = origin_header
        mrc.close()

def main(config_dict):
    cord_data = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_coordinates_ca.txt"
    cord_probs_aa = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_amino_atom_common_emi.txt"
    cord_probs_sec = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_probabilities_sec_emi.txt"
    save_cords = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_cluster_transition_ca.txt"
    save_probs_aa = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_cluster_emission_aa_ca.txt"
    save_probs_sec = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_cluster_emission_sec_ca.txt"
    org_mrc = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/emd_normalized_map.mrc"
    save_cords_mrc = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_cluster_ca.mrc"
    save_cords_pdb = f"{config_dict['test_data_dir']}/{config_dict['density_map_name']}/{config_dict['density_map_name']}_cluster_ca.pdb"

    if os.path.exists(save_cords):
        os.remove(save_cords)

    if os.path.exists(save_probs_aa):
        os.remove(save_probs_aa)

    if os.path.exists(save_probs_sec):
        os.remove(save_probs_sec)

    if os.path.exists(save_cords_pdb):
        os.remove(save_cords_pdb)
        
    print("Running Cluster for: ", cord_data, " with threshold:", config_dict['clustering_threshold'])
    proc_probabilities_aa(cord_probs_aa)
    centroid(cord_data, save_cords, save_probs_aa, save_probs_sec, config_dict['clustering_threshold'])
    generate_pdb_mrc(save_cords=save_cords, save_cords_pdb=save_cords_pdb, save_cords_mrc=save_cords_mrc, org_mrc=org_mrc)
    print("Cluster PDB Done and saved in : ", save_cords_pdb)
    print("Cluster MRC Done and saved in : ", save_cords_mrc)
    return save_cords, save_probs_aa
