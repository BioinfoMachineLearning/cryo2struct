"""
Created on 6 March 2023 01:15 PM
@author: nabin


This script takes in predicted probability file, process it and extracts only ca from it, then finally saves to mrc file.

"""
import mrcfile
from tqdm import tqdm
import numpy as np
import ast
import math
import os

def get_index(cord, origin, voxel):
    return math.ceil(math.floor(cord - abs(origin)) / voxel)


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
        fi.write('A')  # todo: need to change chain id accordingly
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


def extract_ca_from_atom_amino_sec_common_only(pred_atom, pred_amino, pred_sec ,outfilename, outfilename_pdb, density_shape, density_voxel, density_origin, origin):
    # get common from atom, sec and amino pred
    data = np.zeros(density_shape, dtype=np.float32)
    atom_idx = dict()
    atom_sec_idx = dict()
    count = 0
    key_err = 0
    idx_err = 0
    idx_no_err = 0
    x_origin = density_origin[0]
    y_origin = density_origin[1]
    z_origin = density_origin[2]
    x_voxel = density_voxel[0]
    y_voxel = density_voxel[1]
    z_voxel = density_voxel[2]

    with open(pred_atom, 'r') as atom_prob:
        for line in atom_prob:
            line_a = ast.literal_eval(line)
            ca_prob = line_a[2]
            ca_cords = line_a[0]
            x, y, z = ca_cords[0], ca_cords[1], ca_cords[2]
            iz = int(get_index(z, z_origin, z_voxel))
            jy = int(get_index(y, y_origin, y_voxel))
            kx = int(get_index(x, x_origin, x_voxel))
            atom_idx[iz,jy,kx] = ca_prob


    with open(pred_sec, 'r') as sec_prob:
        for line in sec_prob:
            line_a = ast.literal_eval(line)
            ca_prob = 1 - line_a[1]
            ca_cords = line_a[0]
            x, y, z = ca_cords[0], ca_cords[1], ca_cords[2]
            iz = int(get_index(z, z_origin, z_voxel))
            jy = int(get_index(y, y_origin, y_voxel))
            kx = int(get_index(x, x_origin, x_voxel))
            try:
                atom_idx[iz,jy,kx] = np.sqrt(atom_idx[iz,jy,kx] * ca_prob)
                atom_sec_idx[iz,jy,kx] = atom_idx[iz,jy,kx]
            except KeyError:
                key_err += 1


    with open(pred_amino, 'r') as amino_prob:
        for line in amino_prob:
            line_a = ast.literal_eval(line)
            ca_prob = 1 - line_a[1]
            ca_cords = line_a[0]
            x, y, z = ca_cords[0], ca_cords[1], ca_cords[2]
            iz = int(get_index(z, z_origin, z_voxel))
            jy = int(get_index(y, y_origin, y_voxel))
            kx = int(get_index(x, x_origin, x_voxel))
            try:
                atom_sec_idx[iz,jy,kx] = np.sqrt(atom_sec_idx[iz,jy,kx] * ca_prob)
                data[iz,jy,kx] = atom_sec_idx[iz,jy,kx] 
                save(x=x, y=y, z=z, count=count, save_path=outfilename_pdb)
                count += 1
                idx_no_err += 1
            except KeyError:
                key_err += 1
                idx_err += 1
            
    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = origin
        mrc.close()


    print("####################################################")
    print("Atom_Sec_Amino index error", idx_err)
    print("Atom_Sec_Amino NO index error", idx_no_err)
    print("Number of common carbon alphas", idx_no_err)

def extract_ca_from_atom_amino_common_only(pred_atom, pred_amino, outfilename, outfilename_pdb, density_shape, density_voxel, density_origin, origin):
    data = np.zeros(density_shape, dtype=np.float32)
    atom_idx = dict()
    count = 0
    key_err = 0
    idx_err = 0
    idx_no_err = 0
    x_origin = density_origin[0]
    y_origin = density_origin[1]
    z_origin = density_origin[2]
    x_voxel = density_voxel[0]
    y_voxel = density_voxel[1]
    z_voxel = density_voxel[2]

    with open(pred_atom, 'r') as atom_prob:
        for line in atom_prob:
            line_a = ast.literal_eval(line)
            ca_prob = line_a[2]
            ca_cords = line_a[0]
            x, y, z = ca_cords[0], ca_cords[1], ca_cords[2]
            iz = int(get_index(z, z_origin, z_voxel))
            jy = int(get_index(y, y_origin, y_voxel))
            kx = int(get_index(x, x_origin, x_voxel))
            atom_idx[iz,jy,kx] = ca_prob

    with open(pred_amino, 'r') as amino_prob:
        for line in amino_prob:
            line_a = ast.literal_eval(line)
            ca_prob = 1 - line_a[1]
            ca_cords = line_a[0]
            x, y, z = ca_cords[0], ca_cords[1], ca_cords[2]
            iz = int(get_index(z, z_origin, z_voxel))
            jy = int(get_index(y, y_origin, y_voxel))
            kx = int(get_index(x, x_origin, x_voxel))
            try:
                # atom_idx[iz,jy,kx] = np.sqrt(atom_idx[iz,jy,kx] * ca_prob)
                a = np.sqrt(atom_idx[iz,jy,kx] * ca_prob)
                data[iz,jy,kx] = atom_idx[iz,jy,kx] 
                save(x=x, y=y, z=z, count=count, save_path=outfilename_pdb)
                count += 1
                idx_no_err += 1
            except KeyError:
                key_err += 1
                idx_err += 1
            

    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = origin
        mrc.close()

    print("####################################################")
    print("Atom_Amino index error", idx_err)
    print("Atom_Amino NO index error", idx_no_err)
    print("Number of common carbon alphas", idx_no_err)


def extract_ca_from_amino(pred_amino, outfilename, outfilename_pdb, density_shape, density_voxel, density_origin, origin):
    # use this for data visualization, only ca from amino prediction
    data = np.zeros(density_shape, dtype=np.float32)
    count = 0
    idx_err = 0
    idx_no_err = 0
    x_origin = density_origin[0]
    y_origin = density_origin[1]
    z_origin = density_origin[2]
    x_voxel = density_voxel[0]
    y_voxel = density_voxel[1]
    z_voxel = density_voxel[2]
    with open(pred_amino, 'r') as atom_prob:
        for line in atom_prob:
            line_a = ast.literal_eval(line)
            ca_prob = 1 - line_a[1]
            ca_cords = line_a[0]
            x, y, z = ca_cords[0], ca_cords[1], ca_cords[2]
            save(x=x, y=y, z=z, count=count, save_path=outfilename_pdb)
            count += 1
            iz = int(get_index(z, z_origin, z_voxel))
            jy = int(get_index(y, y_origin, y_voxel))
            kx = int(get_index(x, x_origin, x_voxel))
            try:
                data[iz,jy,kx] = ca_prob  
                idx_no_err += 1  
            except:
                idx_err += 1

    print("####################################################")
    print("Amino index error", idx_err)
    print("Amino NO index error", idx_no_err)


    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = origin
        mrc.close()



def extract_ca_from_atom(pred_atom, outfilename, outfilename_pdb, density_shape, density_voxel, density_origin, origin):
    # use this for data visualization, only ca from atom prediction
    data = np.zeros(density_shape, dtype=np.float32)
    idx_err = 0
    idx_no_err = 0
    count = 0
    x_origin = density_origin[0]
    y_origin = density_origin[1]
    z_origin = density_origin[2]
    x_voxel = density_voxel[0]
    y_voxel = density_voxel[1]
    z_voxel = density_voxel[2]
    with open(pred_atom, 'r') as atom_prob:
        for line in atom_prob:
            line_a = ast.literal_eval(line)
            ca_prob = line_a[2]
            ca_cords = line_a[0]
            x, y, z = ca_cords[0], ca_cords[1], ca_cords[2]
            save(x=x, y=y, z=z, count=count, save_path=outfilename_pdb)
            count += 1
            iz = int(get_index(z, z_origin, z_voxel))
            jy = int(get_index(y, y_origin, y_voxel))
            kx = int(get_index(x, x_origin, x_voxel))
            try:
                data[iz,jy,kx] = ca_prob  
                idx_no_err += 1  
            except:
                idx_err += 1
    
    print("####################################################")
    print("Atom index error", idx_err)
    print("Atom NO index error", idx_no_err)
    
    with mrcfile.new(outfilename, overwrite=True) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = x_voxel
        mrc.header.origin = origin
        mrc.close()