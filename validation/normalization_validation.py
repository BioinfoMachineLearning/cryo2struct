
"""
@author: nabin

- Validates the normalization of density maps.
"""

import numpy as np
import os
import mrcfile
from copy import deepcopy


if __name__ == "__main__":
    input_directory = "cryo2struct/root/EMD_0"
    count_error = 0
    count_no_error = 0
    density_maps = os.listdir(input_directory)
    for maps in density_maps:
        normalized_map = f"{input_directory}/{maps}/emd_normalized_map.mrc"
        clean_map = mrcfile.open(normalized_map, mode='r')
        map_data = deepcopy(clean_map.data)
        all_values_within_range = np.all((map_data >= 0) & (map_data <= 1))
        if not all_values_within_range:
            print("Not all values are within the range [ 0, 1].", maps)
            count_error += 1
        else:
            print("Done", maps)
            count_no_error += 1
    print(f"Total maps: {len(density_maps)}. No Error maps : {count_no_error}. Error maps : {count_error}. {input_directory}")
 