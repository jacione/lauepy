import numpy as np
from scipy.spatial.transform import Rotation


def rot(data, axis, angle):
    for n in range(len(axis)):
        r = Rotation.from_euler(axis[n], angle[n], degrees=True).as_matrix()
        data = r @ data
        data = data.T
    return data


def transform(abc_star, phi, chi, theta):
    M_rot = rot(abc_star, ['y', 'z', 'x'], [-theta, chi - 90, -phi])
    out_plane = M_rot @ np.array([0, 1, 0]).T
    in_plane = M_rot @ np.array([1, 0, 0]).T

    return in_plane.T, out_plane.T


def spec_to_sample(hkl, phi, chi, theta):
    new_hkl = rot(hkl.T, ['x', 'z', 'y'], [-phi, -(chi - 90), -theta])

    return new_hkl


def sample_to_spec(hkl, phi, chi, theta):
    new_hkl = rot(hkl.T, ['y', 'z', 'x'], [theta, chi - 90, phi])

    return new_hkl
