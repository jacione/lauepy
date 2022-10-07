import json
import subprocess as sub
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from src.lauepy.disorientation import rmat_2_quat, calc_disorient


def find_possible_twins(working_dir=None, twin_tolerance=None, verbose=True, **kwargs):
    with open(Path(working_dir) / "grains/grains.json") as f:
        grains = json.load(f)

    rot_mats = [np.array(grains[grain]['Rot_mat']) for grain in grains]
    keys = list(grains)
    indices = np.arange(len(keys))
    comb = list(combinations(indices, 2))
    # print(comb)
    q1 = rmat_2_quat([np.array(grains[keys[c[0]]]['Rot_mat']) for c in comb])
    q2 = rmat_2_quat([np.array(grains[keys[c[1]]]['Rot_mat']) for c in comb])
    misorientation = calc_disorient(q1, q2)
    good_ind = np.where(np.absolute(misorientation - 60) < twin_tolerance)[0]
    inds = np.array([comb[g] for g in good_ind])
    if verbose:
        print("Possible twins:")
        print(inds-1)
        print()
    with open(Path(working_dir) / "twins/possible_twins.txt", 'w') as file:
        for i in inds:
            file.write("%d %d \n" % tuple(i))
    return bool(len(inds))


def find_twins(working_dir=None, twin_tolerance=None, **kwargs):
    with open(f"{working_dir}/grains/grains.json", 'r') as f:
        grains = json.load(f)

    keys = list(grains)

    inds = np.array(np.loadtxt(f"{working_dir}/twins/possible_twins.txt", dtype=np.int32))
    if np.shape(inds) == (2,):
        inds = [inds]
    e1 = np.array(
        [Rotation.from_matrix(np.array(grains[keys[i[0]]]['Rot_mat']).T).as_euler('ZXZ', degrees=True) for i in inds])
    e2 = np.array(
        [Rotation.from_matrix(np.array(grains[keys[i[1]]]['Rot_mat']).T).as_euler('ZXZ', degrees=True) for i in inds])
    es = np.hstack([e1, e2])
    id_pairs = [(keys[i[0]], keys[i[1]]) for i in inds]
    angle_path = f"{working_dir}/twins/angle_list.txt"
    with open(angle_path, 'w') as file:
        file.writelines('angles')
        for e in es:
            file.write('\n%s %s %s %s %s %s' % tuple(e))

    angle_calculator = Path(__file__).resolve().parents[0] / 'a.out'
    sub.run([f'{angle_calculator}', '0', '1', '1', '0', angle_path], stdout=sub.DEVNULL, stderr=sub.DEVNULL)
    data = np.loadtxt('min-misor-angles.txt', skiprows=1)
    if len(list(data.shape)) == 1:
        data = np.array([data])
    data = data[:, 2:6]
    #     data = np.concatenate(data,axis=1)
    s1 = [np.array(grains[keys[i[0]]]['Spec_Orientation']).ravel() for i in inds]
    s2 = [np.array(grains[keys[i[1]]]['Spec_Orientation']).ravel() for i in inds]
    ss = np.hstack([s1, s2])
    twin_orientations = data

    twin_orientations = [t for t in twin_orientations
                         if np.sum(np.absolute(np.absolute(t[-3:]) - 0.577)) < twin_tolerance]
    for twin in twin_orientations:
        twin[-3:] = twin[-3:] / twin[-3:].min()
    with open(f"{working_dir}/twins/hiconf_twins.txt", 'w') as file:
        file.writelines('ID_A ID_B Spec_Ori1 Spec_Ori2 mis ax1 ax2 ax3')
        for i, tw in enumerate(twin_orientations):
            file.write('\n%s %s %.2f %.3f %.3f %.3f' % (id_pairs[i] + tuple(tw)))
    return


def cleanup_directory():
    intermediate_files = [
        "min-misor-angles.txt",
        "misor-angles-pairs.txt",
        "OutputForAnalysis.txt",
        "rexgbs_Eulers.wts",
        "rexgbs_misors.mdf",
        "rexgbs_misors_forward.mdf",
        "rexgbs_output.txt"
    ]
    main_dir = Path(__file__)
    for f in intermediate_files:
        (main_dir.parents[1] / f"lauepy_scripts/{f}").unlink()
