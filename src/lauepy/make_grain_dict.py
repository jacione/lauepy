import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
from sklearn import cluster

from src.lauepy.disorientation import calc_disorient, rmat_2_quat
from src.lauepy.write_macro import grain_to_macro


def cluster_grains(config):
    print('Grouping patterns into grains...')

    working_dir = config['working_dir']
    for p in Path(f"{working_dir}/grains").iterdir():
        p.unlink()
    for p in Path(f"{working_dir}/macros").iterdir():
        p.unlink()
    with open(f'{working_dir}/peaks/patterns.json', 'r') as f:
        pattern_dict = json.load(f)

    patterns = [p for p in pattern_dict.values()]

    print(rmat_2_quat([p["Rot_mat"] for p in patterns]))

    orients = np.array([p["Spec_Orientation"] for p in patterns])
    orients = np.sort(np.abs(orients), axis=-1)
    orients /= orients[:, :, -1].reshape((*orients[:, :, -1].shape, 1))
    orients = orients.reshape(-1, 6)  # The last element of each normalized hkl vector is constrained to unity.

    clustering = cluster.DBSCAN(eps=0.0005, min_samples=1, metric="cosine")
    # clustering = cluster.AgglomerativeClustering(n_clusters=None, distance_threshold=0.02, linkage='ward')
    clustering.fit(orients)
    labels = clustering.labels_
    print(f"Found {len(np.unique(labels))} grains")
    print(f"Ignored {np.sum(labels == -1)} patterns")

    grains = {f"grain_{i}": None for i in np.unique(labels) if i != -1}
    for pattern, label in zip(patterns, labels):
        if label == -1:
            continue
        if pattern["Center_Frame"] == -1:
            grains["substrate"] = append_to_grain(pattern, None)
            del grains[f"grain_{label}"]
        else:
            grains[f"grain_{label}"] = append_to_grain(pattern, grains[f"grain_{label}"])

    for grain in grains.values():
        grain["orientation"] = np.average(grain["orientation"], axis=0, weights=grain["num_peaks"])
        grain["num_frames"] = len(grain["num_frames"])
        grain["position"] = np.average(grain["position"], axis=0, weights=grain["num_peaks"])
        grain["rms_error"] = np.mean(grain["rms_error"])
        grain["num_peaks"] = np.mean(grain["num_peaks"])
    print_clusters(grains, config)
    map_clusters(grains, config)


def append_to_grain(pattern, grain):
    if grain is None:
        grain = {
            "orientation": [np.sort(np.abs(pattern["Spec_Orientation"]), axis=-1).tolist()],
            "num_frames": [pattern["Center_Frame"]],
            "num_peaks": [pattern["Num_Peaks"]],
            "position": [pattern["Pos"]],
            "rms_error": [pattern["RMS"]]
        }
    else:
        grain["orientation"].append(np.sort(np.abs(pattern["Spec_Orientation"]), axis=-1).tolist())
        grain["num_frames"].append(pattern["Center_Frame"])
        grain["num_peaks"].append(pattern["Num_Peaks"])
        grain["position"].append(pattern["Pos"])
        grain["rms_error"].append(pattern["RMS"])
    return grain


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

    num_grains = len([key for key in grains.keys() if key.startswith("grain")])
    num_frames = [len(grains[f"grain_{num}"]["Frames"]) for num in range(1, num_grains + 1)]
    sorted_grains = {key: val for key, val in grains.items() if not key.startswith("grain")}
    for i, j in enumerate(np.flip(np.argsort(num_frames))):
        sorted_grains[f"grain_{i}"] = grains[f"grain_{j + 1}"]

    with open(f'{working_dir}/grains/grains.json', 'w') as json_file:
        json.dump(sorted_grains, json_file)

    grain_to_macro(config)
    with open(f'{working_dir}/peaks/patterns.json', 'w') as json_file:
        json.dump(pattern_dict, json_file)
    print_grains(sorted_grains, config)
    map_grains(config)
    return


def print_grains(grains, config):
    print('################### GRAIN DICT #######################')
    print(f'{"Grain":>10}{"RMS":>8}{"Peaks":>7}{"Frames":>8}   {"Position":<16}{"HKL-in":<16}{"HKL-out":<16}')
    for grain in grains:
        g = grains[grain]
        num_frames = len(g['Frames'])
        if num_frames > config["grain_threshold"]:
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


def print_clusters(grains, config):
    s = '################### GRAINS #######################\n'
    s += f'{"Grain":>10}{"RMS":>8}{"Peaks":>7}{"Frames":>8}   {"Position":<16}{"HKL-in":<16}{"HKL-out":<16}\n'
    for grain in grains:
        g = grains[grain]
        # if g['num_frames'] > config["grain_threshold"]:
        #     continue
        hkl_in = [round(x) for x in g['orientation'][0]]
        hkl_out = [round(x) for x in g['orientation'][1]]
        s += f"{grain:>10}"
        s += f"{g['rms_error']:8.3f}"
        s += f"{g['num_peaks']:7.1f}"
        s += f"{g['num_frames']:>8}   "
        s += f"[{g['position'][0]:6.2f},{g['position'][1]:6.2f}] "
        s += f"[{hkl_in[0]:>3}, {hkl_in[1]:>3}, {hkl_in[2]:>3}] "
        s += f"[{hkl_out[0]:>3}, {hkl_out[1]:>3}, {hkl_out[2]:>3}]\n"
    (Path(config['working_dir'])/'grains/grain_output.txt').write_text(s)
    if config['verbose']:
        print(s)


def map_clusters(grains, config):
    working_dir = config["working_dir"]

    # Load and select the grains to show
    # with open(f'{working_dir}/grains/grains.json', 'r') as f:
    #     grains = {
    #         grain[6:]: g["COM"]
    #         for grain, g in json.load(f).items()
    #         if (1 < len(g['Frames']) or g['Avg_Peaks'] > 3.01)
    #            and len(g["Frames"]) < config["grain_threshold"]
    #            and grain != "substrate"
    #     }
    grains = {
        grain[6:]: g["position"]
        for grain, g in grains.items()
        if (1 < g['num_frames'] or g['num_peaks'] > 3.01)
           and g["num_frames"] < config["grain_threshold"]
           and grain != "substrate"
    }
    coords = np.array([c for c in grains.values()])

    # Transform the labx/labz coordinates into the sample frame of reference (top-down, beam travelling up)
    with open(f'{working_dir}/peaks/peaks.json', 'r') as f:
        phi, chi, theta = json.load(f)["info"]["angles"]
    alpha = np.abs(np.deg2rad(phi))
    beta = np.abs(np.deg2rad(chi - 90))
    theta = np.abs(np.deg2rad(theta))
    grazing = np.pi / 2 + alpha * np.cos(theta) + beta * np.sin(theta)
    coords[:, 1] *= np.tan(grazing)

    # Scatterplot the COM grain positions
    plt.figure()
    plt.xlabel('microns')
    plt.ylabel('microns')
    plt.gca().set_aspect('equal')
    plt.scatter(*coords.T, c='b')
    ts = []
    for num, c in zip(grains, coords):
        ts.append(plt.text(*c, num, size=15))
    adjust_text(ts, x=coords[:, 0], y=coords[:, 1], force_points=0.25, lim=50)
    plt.savefig(f"{working_dir}/grains/scatterplot.png", dpi=300)
    if config["show_plots"]:
        plt.show()


def map_grains(config):
    working_dir = config["working_dir"]

    # Load and select the grains to show
    with open(f'{working_dir}/grains/grains.json', 'r') as f:
        grains = {
            grain[6:]: g["COM"]
            for grain, g in json.load(f).items()
            if (1 < len(g['Frames']) or g['Avg_Peaks'] > 3.01)
               and len(g["Frames"]) < config["grain_threshold"]
               and grain != "substrate"
        }
    coords = np.array([c for c in grains.values()])

    # Transform the labx/labz coordinates into the sample frame of reference (top-down, beam travelling up)
    with open(f'{working_dir}/peaks/peaks.json', 'r') as f:
        phi, chi, theta = json.load(f)["info"]["angles"]
    alpha = np.abs(np.deg2rad(phi))
    beta = np.abs(np.deg2rad(chi - 90))
    theta = np.abs(np.deg2rad(theta))
    grazing = np.pi / 2 + alpha * np.cos(theta) + beta * np.sin(theta)
    coords[:, 1] *= np.tan(grazing)

    # Scatterplot the COM grain positions
    plt.figure()
    plt.xlabel('microns')
    plt.ylabel('microns')
    plt.gca().set_aspect('equal')
    plt.scatter(*coords.T, c='b')
    ts = []
    for num, c in zip(grains, coords):
        ts.append(plt.text(*c, num, size=15))
    adjust_text(ts, x=coords[:, 0], y=coords[:, 1], force_points=0.25)
    plt.savefig(f"{working_dir}/grains/scatterplot.png", dpi=300)
    if config["show_plots"]:
        plt.show()
