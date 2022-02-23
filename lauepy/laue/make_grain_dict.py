import json

import numpy as np

import lauepy.laue.overlay_peaks as op
from lauepy.laue.disorientation import calc_disorient, rmat_2_quat
from lauepy.laue.write_specorient import write_orient


def make_grain_dict(config):
    working_dir = config['working_dir']
    with open(f'{working_dir}/peaks/patterns.json', 'r') as f:
        pattern_dict = json.load(f)
    grains = {}
    count = 1

    for patt in pattern_dict:
        if pattern_dict[patt]['Center_Frame'] == -1:
            grains['substrate'] = {
                'Rot_mat': pattern_dict[patt]['Rot_mat'],
                'Spec_Orientation': pattern_dict[patt]['Spec_Orientation'],
                'Patts': [patt],
                'Frames': [pattern_dict[patt]['Center_Frame']],
                'Positions': [pattern_dict[patt]['Pos']],
                'RMS': [float(np.round(pattern_dict[patt]['RMS'], 2))],
                'Dist': [float(np.round(pattern_dict[patt]['Dist'], 3))],
                'Num_Peaks': [float(np.round(pattern_dict[patt]['Num_Peaks'], 1))],
                'Patterns': [patt],
                'Count': [pattern_dict[patt]['Count']]
            }
            continue
        match = False
        if len(list(grains)) > 0:
            pat_rot = [np.array(pattern_dict[patt]['Rot_mat']) for _ in grains]
            grain_rot = [np.array(grains[grain]['Rot_mat']) for grain in grains]
            misorientations = calc_disorient(rmat_2_quat(pat_rot), rmat_2_quat(grain_rot))
            for key, grain in enumerate(grains):
                if misorientations[key] < config['grain_tolerance']:
                    if pattern_dict[patt]['Center_Frame'] not in grains[grain]['Frames']:
                        grains[grain]['Frames'].append(pattern_dict[patt]['Center_Frame'])
                        grains[grain]['Positions'].append(pattern_dict[patt]['Pos'])
                    grains[grain]['Patts'].append(patt)
                    grains[grain]['RMS'].append(float(np.round(pattern_dict[patt]['RMS'], 2)))
                    grains[grain]['Dist'].append(float(np.round(pattern_dict[patt]['Dist'], 3)))
                    grains[grain]['Num_Peaks'].append(float(np.round(pattern_dict[patt]['Num_Peaks'], 1)))
                    grains[grain]['Patterns'].append(patt)
                    grains[grain]['Count'].append(pattern_dict[patt]['Count'])
                    match = True
        if not match:
            grains[f'grain_{count}'] = {
                'Rot_mat': pattern_dict[patt]['Rot_mat'],
                'Spec_Orientation': pattern_dict[patt]['Spec_Orientation'],
                'Patts': [patt],
                'Frames': [pattern_dict[patt]['Center_Frame']],
                'Positions': [pattern_dict[patt]['Pos']],
                'RMS': [float(np.round(pattern_dict[patt]['RMS'], 2))],
                'Dist': [float(np.round(pattern_dict[patt]['Dist'], 3))],
                'Num_Peaks': [float(np.round(pattern_dict[patt]['Num_Peaks'], 1))],
                'Patterns': [patt],
                'Count': [pattern_dict[patt]['Count']]
            }
            count += 1
    for grain in grains:
        grains[grain]['Avg_RMS'] = float(np.round(np.mean(grains[grain]['RMS']), 2))
        grains[grain]['Avg_Dist'] = float(np.round(np.mean(grains[grain]['Dist']), 4))
        grains[grain]['Avg_Peaks'] = float(np.round(np.mean(grains[grain]['Num_Peaks']), 1))
        grains[grain]['Avg_Count'] = float(np.round(np.mean(grains[grain]['Count']), 1))
        #         an_array[numpy.argsort(an_array[:, 1])]
        pos = np.array(grains[grain]['Positions'])
        grains[grain]['COM'] = pos[np.argsort(pos[:, 0])][len(pos) // 2].tolist()

    with open(f'{working_dir}/grains/grains.json', 'w') as json_file:
        json.dump(grains, json_file)
    # if not os.path.exists(working_dir):
    #     os.mkdir(working_dir)
    # else:
    #     shutil.rmtree(working_dir)
    #     os.mkdir(working_dir)
    # if os.path.exists(working_dir + '/seg_stack.tif'):
    #     os.remove(working_dir + '/seg_stack.tif')
    # copyfile('seg_stack.tif', working_dir + '/seg_stack.tif')

    for key, grain in grains.items():
        if key == 'substrate':
            continue
        for patt in grain['Patts']:
            pattern_dict[patt]['Grain'] = grain
        # print("g['Patts']",g['Patts'],"g['Positions']",g['Positions'],"g['COM']",g['COM'])
        # med_patt = g['Patts'][g['Positions'].index(g['COM'])]
        # print("med_patt",med_patt,"grain",grain)
        draw_patts = grain['Patts']
        for pattern in draw_patts:
            op.overlay(config, pattern_dict[pattern], key)
            # print(output_directory,len(grains))
    write_orient(config)
    with open(f'{working_dir}/peaks/patterns.json', 'w') as json_file:
        json.dump(pattern_dict, json_file)
    return
