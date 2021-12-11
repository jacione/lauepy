import json
import os
import shutil
from shutil import copyfile

import numpy as np

import lauepy.laue.overlay_peaks as op
from lauepy.laue.disorientation import calc_disorient, rmat_2_quat
from lauepy.laue.write_specorient import write_orient as wo


def make_grain_dict(output_directory, pattern_dict_file='pattern_dict.json', grain_dict_file='grain_dict.json',
                    cryst_path='crystal_params.json', mis_tol=0.5, plot=True, threshold=1e4,
                    det_params='det_params.json', frame_stepsize=1):
    with open(pattern_dict_file) as f:
        pattern_dict = json.load(f)
    with open(grain_dict_file) as f:
        grains = json.load(f)
    count = 1

    for patt in pattern_dict:

        match = False
        if len(list(grains)) > 0:
            pat_rot = [np.array(pattern_dict[patt]['Rot_mat']) for _ in grains]
            grain_rot = [np.array(grains[grain]['Rot_mat']) for grain in grains]
            misorientations = calc_disorient(rmat_2_quat(pat_rot), rmat_2_quat(grain_rot))
        else:
            continue

        for i, grain in enumerate(grains):

            if misorientations[i] < mis_tol:
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
            grains['grain_%s' % count] = {}
            grains['grain_%s' % count]['Rot_mat'] = pattern_dict[patt]['Rot_mat']
            grains['grain_%s' % count]['Spec_Orientation'] = pattern_dict[patt]['Spec_Orientation']
            grains['grain_%s' % count]['Patts'] = [patt]
            grains['grain_%s' % count]['Frames'] = [pattern_dict[patt]['Center_Frame']]
            grains['grain_%s' % count]['Positions'] = [pattern_dict[patt]['Pos']]
            grains['grain_%s' % count]['RMS'] = [float(np.round(pattern_dict[patt]['RMS'], 2))]
            grains['grain_%s' % count]['Dist'] = [float(np.round(pattern_dict[patt]['Dist'], 3))]
            grains['grain_%s' % count]['Num_Peaks'] = [float(np.round(pattern_dict[patt]['Num_Peaks'], 1))]
            grains['grain_%s' % count]['Patterns'] = [patt]
            grains['grain_%s' % count]['Count'] = [pattern_dict[patt]['Count']]
            count += 1
    for grain in grains:
        grains[grain]['Avg_RMS'] = float(np.round(np.mean(grains[grain]['RMS']), 2))
        grains[grain]['Avg_Dist'] = float(np.round(np.mean(grains[grain]['Dist']), 4))
        grains[grain]['Avg_Peaks'] = float(np.round(np.mean(grains[grain]['Num_Peaks']), 1))
        grains[grain]['Avg_Count'] = float(np.round(np.mean(grains[grain]['Count']), 1))
        #         an_array[numpy.argsort(an_array[:, 1])]
        pos = np.array(grains[grain]['Positions'])
        grains[grain]['COM'] = pos[np.argsort(pos[:, 0])][len(pos) // 2].tolist()

    with open(grain_dict_file, 'w') as json_file:
        json.dump(grains, json_file)
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    else:
        shutil.rmtree(output_directory)
        os.mkdir(output_directory)
    if os.path.exists(output_directory + '/seg_stack.tif'):
        os.remove(output_directory + '/seg_stack.tif')
    copyfile('seg_stack.tif', output_directory + '/seg_stack.tif')

    for grain in grains:

        g = grains[grain]
        for patt in g['Patts']:
            pattern_dict[patt]['Grain'] = grain
        # print("g['Patts']",g['Patts'],"g['Positions']",g['Positions'],"g['COM']",g['COM'])
        # med_patt = g['Patts'][g['Positions'].index(g['COM'])]
        # print("med_patt",med_patt,"grain",grain)
        draw_patts = g['Patts']
        for draw_patt in draw_patts:
            op.overlay(output_directory, draw_patt, grain, threshold=threshold, reduction='reduction.yml',
                       pat_dict='pattern_dict.json', plot=plot, det_params=det_params, cryst_params=cryst_path,
                       frame_step=frame_stepsize)
            # print(output_directory,len(grains))
    wo(output_directory, cryst_path=cryst_path, grain_path=grain_dict_file)
    with open(pattern_dict_file, 'w') as json_file:
        json.dump(pattern_dict, json_file)
    return
