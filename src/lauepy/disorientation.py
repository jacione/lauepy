from scipy.spatial.transform import Rotation

import numpy as np
import cupy as cp


# convert euler angles to quaternions
def rmat_2_quat(rmat):
    r = Rotation.from_matrix(rmat)

    # Invert to match the massif convention, where the Bunge Euler/Rotation is a
    # transformation from sample to crystal frame.
    r1 = r.inv()
    quat = r1.as_quat()

    for num, val in enumerate(quat):
        if val[3] < 0:  # q1,q2,q3,q0 format
            quat[num] = -1.0 * quat[num]
    quat = np.roll(quat, 1, axis=1)

    rmat = r1.as_matrix()
    rmat_inv = r.as_matrix()

    return quat


# compute quaternion product
def QuadProd(p, q):
    p = cp.array(p)
    q = cp.array(q)
    p0 = cp.reshape(p[:, 0], (p[:, 0].shape[0], 1))
    q0 = cp.reshape(q[:, 0], (q[:, 0].shape[0], 1))

    l = cp.sum(p[:, 1:]*q[:, 1:], 1)

    prod1 = (p0*q0)[:, 0] - l

    prod2 = p0*q[:, 1:] + q0*p[:, 1:] + cp.cross(p[:, 1:], q[:, 1:])
    m = cp.transpose(cp.stack([prod1, prod2[:, 0], prod2[:, 1], prod2[:, 2]]))

    return m


# invert quaternion
def invQuat(p):
    q = cp.transpose(cp.stack([-p[:, 0], p[:, 1], p[:, 2], p[:, 3]]))
    return q


# calculate the disorientation between two sets of quaternions ps and qs
def calc_disorient(y_true, y_pred):
    y_true = cp.array(y_true)
    y_pred = cp.array(y_pred)

    # sort quaternion for cubic symmetry trick
    p = cp.sort(cp.abs(QuadProd(invQuat(y_true), y_pred)))

    # calculate last component of two other options
    p1 = (p[:, 2] + p[:, 3]) / 2 ** (1 / 2)
    p2 = (p[:, 0] + p[:, 1] + p[:, 2] + p[:, 3]) / 2
    vals = cp.transpose(cp.stack([p[:, -1], p1, p2]))
    # pick largest value and find angle

    max_val = cp.amax(vals, axis=1)
    mis = (2 * cp.arccos(max_val))

    return cp.degrees(replacenan(mis)).get()


# replace NaNs with zeros
def replacenan(t):
    t[cp.isnan(t)] = 0
    return t
