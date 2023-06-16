
"""
@author: nabin

- Runs ChimeraX in no-GUI mode to resample map.
"""
import sys
import subprocess
import os


def execute(input_path, chimera_path):
    """
    creates resampling script and executes them in Chimera
    :return:
    """
    map_names = [fn for fn in os.listdir(input_path) if not fn.endswith(".ent") if not fn.endswith(".DS_Store")]
    if map_names is None:
        print("### Please check the directory!! No input files present in", input_path, '###')
        exit()
    for maps in range(len(map_names)):
        path = os.path.join(input_path, map_names[maps])
        emd_map = [e for e in os.listdir(path) if e.endswith(".map")]
        for density_maps in range(len(emd_map)):
            chimera_scripts = open('resample.cxc', 'w')
            chimera_scripts.write('open ' + path + '/' + emd_map[density_maps] + '\n'
                                  'vol resample #1 spacing 1.0 \n'
                                  'save ' + path + '/' + emd_map[density_maps].split("_")[0] + '_resampled_map.mrc' + ' model #2 \n'
                                  'exit')

            chimera_scripts.close()
            script_finished = False
            while not script_finished:
                try:
                    subprocess.run([chimera_path, '--nogui', chimera_scripts.name])
                    script_finished = True
                except FileNotFoundError as error:
                    raise error
            os.remove(chimera_scripts.name)
            print(f'### Resampled {emd_map[density_maps]} and saved on new grid with voxel size of {1} ###')


if __name__ == "__main__":
    input_path = sys.argv[1]
    # output_path = sys.argv[2]

    if len(sys.argv) > 2:
        chimera_path = sys.argv[2]
    else:
        chimera_path = '/usr/bin/chimerax'

    execute(input_path, chimera_path)
    print("Resampling Complete!")

