
"""
@author: nabin

- Validates the voxel size of density maps.
"""

import mrcfile
import numpy as np
import os

count_success = 0
count_failed = 0

def check_voxel_size(density_map):

    global count_success
    global count_failed

    density_map_mrc = mrcfile.open(density_map, mode='r')
    # Get voxel dimensions from MRC file (unit in Angstrom)
    x_voxel = density_map_mrc.voxel_size['x']
    y_voxel = density_map_mrc.voxel_size['y']
    z_voxel = density_map_mrc.voxel_size['z']

    if x_voxel == 1 and y_voxel == 1 and z_voxel == 1:
        count_success += 1
        print(f"Resampling validation successful for density map {density_map}. Voxel dimensions are 1. ")
    else:
        print(f"Resampling validation failed for density map {density_map}. Voxel dimensions are {x_voxel, y_voxel, z_voxel}")
        count_failed += 1 

if __name__ == "__main__":
    input_directory = "cryo2struct/root/test"
    density_maps = os.listdir(input_directory)
    for maps in density_maps:
        resampled_map = f"{input_directory}/{maps}/emd_resampled_map.mrc"
        check_voxel_size(resampled_map)
    print(f"Total COUNT : {len(density_maps)}")
    print(f"Total SUCCESS : {count_success}")
    print(f"Total FAILED : {count_failed}")
