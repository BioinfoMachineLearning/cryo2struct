"""
@author: nabin

- Normalizes map with 95 percentile, to change the percentile value modify line number 27.
"""
import sys
import mrcfile
from copy import deepcopy
import numpy as np
import os


def execute(inputs):
    count = 0
    map_names = [fn for fn in os.listdir(inputs) if not fn.endswith(".ent")]
    
    for maps in range(len(map_names)):
        if map_names[maps] == "emd_resampled_map.mrc":
            resample_map = map_names[maps]
            os.chdir(inputs)
            print(inputs)
            clean_map = mrcfile.open(resample_map, mode='r')
            map_data = deepcopy(clean_map.data)
            # normalize with percentile value
            print("### Normalizing with 95-percentile for ", resample_map, " ###")
            try:
                percentile = np.percentile(map_data[np.nonzero(map_data)], 95)
                map_data /= percentile
            except IndexError as error:
                count += 1

            # set low valued data to 0
            print("### Setting all values < 0 to 0 for ", resample_map, " ###")
            map_data[map_data < 0] = 0
            print("### Setting all values > 1 to 1 for ", resample_map, " ###")
            map_data[map_data > 1] = 1
            with mrcfile.new(map_names[maps].split("_")[0] + "_normalized_map.mrc", overwrite=True) as mrc:
                mrc.set_data(map_data)
                mrc.voxel_size = 1
                mrc.header.origin = clean_map.header.origin
                mrc.close()
            print("### Wrote file to ", inputs, " ###")
    print("The number of non normalized index: ", count)


if __name__ == "__main__":
    input_path = sys.argv[1]
    maps = [fn for fn in os.listdir(input_path) if not fn.endswith(".DS_Store")]
    for m in range(len(maps)):
        execute(input_path +'/' + maps[m])
    print("Normalization Complete!")
