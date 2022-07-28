import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text

import src.lauepy.overlay_peaks as op
from src.lauepy.disorientation import calc_disorient, rmat_2_quat
from src.lauepy.write_macro import grain_to_macro


def make_grain_dict(config):
    print('Grouping patterns into grains...')

    working_dir = config['working_dir']
    for p in Path(f"{working_dir}/grains").iterdir():
        p.unlink()
    for p in Path(f"{working_dir}/macros").iterdir():
        p.unlink()
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
        pos = np.array(grains[grain]['Positions'])
        grains[grain]['COM'] = np.average(pos, axis=0, weights=np.array(grains[grain]['Num_Peaks'])).tolist()

    with open(f'{working_dir}/grains/grains.json', 'w') as json_file:
        json.dump(grains, json_file)

    for key, grain in grains.items():
        if key == 'substrate':
            continue
        for patt in grain['Patts']:
            pattern_dict[patt]['Grain'] = grain
        draw_patts = grain['Patts']
        for pattern in draw_patts:
            op.overlay(config, pattern_dict[pattern], key)
    grain_to_macro(config)
    with open(f'{working_dir}/peaks/patterns.json', 'w') as json_file:
        json.dump(pattern_dict, json_file)
    if config['verbose']:
        print_grains(grains, config)
    if config['show_plots']:
        map_grains(config)
    return


def print_grains(grains, config):
    print('################### GRAIN DICT #######################')
    print(f'{"Grain":>10}{"RMS":>8}{"Peaks":>7}{"Frames":>8}   {"Position":<16}{"HKL-in":<16}{"HKL-out":<16}')
    for grain in grains:
        g = grains[grain]
        num_frames = len(g['Frames'])
        if (num_frames <= 1 and g['Avg_Peaks'] <= 3.01) or num_frames > config["grain_threshold"]:
            continue
        hkl_in = [round(x) for x in g['Spec_Orientation'][0]]
        hkl_out = [round(x) for x in g['Spec_Orientation'][1]]
        print(f"{grain:>10}"
              f"{g['Avg_RMS']:>8}"
              f"{g['Avg_Peaks']:>7}"
              f"{num_frames:>8}   "
              f"[{np.around(g['COM'][0], 2):>6},{np.around(g['COM'][1], 2):>6}] "
              f"[{hkl_in[0]:>3}, {hkl_in[1]:>3}, {hkl_in[2]:>3}] "
              f"[{hkl_out[0]:>3}, {hkl_out[1]:>3}, {hkl_out[2]:>3}]"
              )


def map_grains(config):
    working_dir = config["working_dir"]

    # Load and select the grains to show
    with open(f'{working_dir}/grains/grains.json', 'r') as f:
        grains = {
            grain[6:]: g["COM"]
            for grain, g in json.load(f).items()
            if (1 < len(g['Frames']) < config["grain_threshold"] or g['Avg_Peaks'] > 3.01) and grain != "substrate"
        }
    coords = np.array([c for c in grains.values()])

    # Transform the labx/labz coordinates into the sample frame of reference (top-down, beam travelling up)
    with open(f'{working_dir}/peaks/peaks.json', 'r') as f:
        phi, chi, theta = json.load(f)["info"]["angles"]
    alpha = np.abs(np.deg2rad(phi))
    beta = np.abs(np.deg2rad(chi-90))
    theta = np.abs(np.deg2rad(theta))
    grazing = np.pi/2 + alpha*np.cos(theta) + beta*np.sin(theta)
    coords[:, 1] *= np.tan(grazing)

    # Scatterplot the grain positions
    plt.figure()
    plt.xlabel('microns')
    plt.ylabel('microns')
    plt.gca().set_aspect('equal')
    plt.scatter(*coords.T, c='b')
    ts = []
    for num, c in zip(grains, coords):
        ts.append(plt.text(*c, num, size=15))
    adjust_text(ts, x=coords[:, 0], y=coords[:, 1], force_points=0.2)
    plt.savefig(f"{working_dir}/grains/scatterplot.png", dpi=300)
    plt.show()
