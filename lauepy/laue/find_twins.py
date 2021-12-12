import json
import subprocess as sub
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from lauepy.laue.disorientation import rmat_2_quat, calc_disorient


def find_possible_twins(grain_dict_path, ang_tol=0.1):
    with open(grain_dict_path) as f:
        grains = json.load(f)

    rot_mats = [np.array(grains[grain]['Rot_mat']) for grain in grains]
    keys = list(grains)
    indices = np.arange(len(keys))
    comb = list(combinations(indices, 2))
    # print(comb)
    q1 = rmat_2_quat([np.array(grains[keys[c[0]]]['Rot_mat']) for c in comb])
    q2 = rmat_2_quat([np.array(grains[keys[c[1]]]['Rot_mat']) for c in comb])
    misorientation = calc_disorient(q1, q2)
    good_ind = np.where(np.absolute(misorientation - 60) < ang_tol)[0]
    inds = np.array([comb[g] for g in good_ind])
    print(len(inds))
    with open('possible_twins.txt', 'w') as file:
        for i in inds:
            file.write("%d %d \n" % tuple(i))
    return


def find_twins(in_dict, in_inds, out_file, deg=0, ax_tol=0.1):
    with open(in_dict) as f:
        grains = json.load(f)

    keys = list(grains)

    inds = np.array(np.loadtxt(in_inds, dtype=np.int32))
    print(inds.shape)
    if np.shape(inds) == (2,):
        inds = [inds]
    e1 = np.array(
        [Rotation.from_matrix(np.array(grains[keys[i[0]]]['Rot_mat']).T).as_euler('ZXZ', degrees=True) for i in inds])
    e2 = np.array(
        [Rotation.from_matrix(np.array(grains[keys[i[1]]]['Rot_mat']).T).as_euler('ZXZ', degrees=True) for i in inds])
    es = np.hstack([e1, e2])
    print(es.shape)
    id_pairs = [(keys[i[0]], keys[i[1]]) for i in inds]
    with open('angle_list.txt', 'w') as file:
        file.writelines('angles')
        for e in es:
            file.write('\n%s %s %s %s %s %s' % tuple(e))

    angle_calculator = Path(__file__).resolve().parents[1] / 'a.out'
    sub.run([f'{angle_calculator}', '0', '1', '1', '0', 'angle_list.txt'])
    data = np.loadtxt('min-misor-angles.txt', skiprows=1)
    if len(list(data.shape)) == 1:
        data = np.array([data])
    data = data[:, 2:6]
    print(data.shape)
    #     data = np.concatenate(data,axis=1)
    s1 = [np.array(grains[keys[i[0]]]['Spec_Orientation']).ravel() for i in inds]
    s2 = [np.array(grains[keys[i[1]]]['Spec_Orientation']).ravel() for i in inds]
    ss = np.hstack([s1, s2])
    print(ss.shape)
    twin_orientations = data

    twin_orientations = [t for t in twin_orientations if np.sum(np.absolute(np.absolute(t[-3:]) - 0.577)) < ax_tol]
    for twin in twin_orientations:
        twin[-3:] = twin[-3:] / twin[-3:].min()
    with open(out_file, 'w') as file:
        file.writelines('ID_A ID_B Spec_Ori1 Spec_Ori2 mis ax1 ax2 ax3')
        for i, tw in enumerate(twin_orientations):
            file.write('\n%s %s %.2f %.3f %.3f %.3f' % (id_pairs[i] + tuple(tw)))
    return
