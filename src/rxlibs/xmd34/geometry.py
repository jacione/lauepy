# -*- coding: utf-8 -*-
"""
Functions & classes for handling detector geometry calculations,
compatible with existing 34ID-E definitions.

Todo:
    * doc of DetectorGeometry

Author : Xu, Ruqing
Created: Jan 2016
"""

import numpy as np

_MIN_RADIAN = 2.0 ** (-26)  # smallest angle in radian that can make numerical difference
# between its cosine (computed in double precision) and 1.0
_COS_MIN_RADIAN = np.cos(_MIN_RADIAN)


def rot_vec_to_matrix(r_vec):
    """
    return rotation matrix based on the Rotation vector,
    formula from http://mathworld.wolfram.com/RodriguesRotationFormula.html
    """
    ang = np.linalg.norm(r_vec)
    if ang > _MIN_RADIAN:
        (wx, wy, wz) = r_vec / ang
    else:
        wx = wy = wz = 0.0
    c = np.cos(ang)
    s = np.sin(ang)
    r_matrix = np.array([[c + (1 - c) * wx ** 2, wx * wy * (1 - c) - wz * s, wy * s + wx * wz * (1 - c)],
                         [wz * s + wx * wy * (1 - c), c + wy * wy * (1 - c), -wx * s + wy * wz * (1 - c)],
                         [-wy * s + wx * wz * (1 - c), wx * s + wy * wz * (1 - c), c + wz * wz * (1 - c)]])

    return r_matrix


def rot_matrix_to_vec(rmatrix):
    """
    return rotation vector from rotation matrix,
    input matrix should be a rotation matrix, no checks performed
    """
    mat = np.array(rmatrix)
    assert mat.shape == (3, 3)
    cosangle = (np.trace(mat) - 1) / 2
    if cosangle > _COS_MIN_RADIAN:
        # 0-degree case,
        return np.array([0, 0, 0], dtype='f8')
    elif cosangle < -_COS_MIN_RADIAN:
        # 180-degree case, assume z positive
        x = np.sqrt(abs(mat[0, 0] + 1) / 2)
        y = np.sqrt(abs(mat[1, 1] + 1) / 2)
        z = np.sqrt(abs(mat[2, 2] + 1) / 2)
        x = -x if (mat[0, 2] + mat[2, 0]) < 0 else x
        y = -y if (mat[1, 2] + mat[2, 1]) < 0 else y
    else:
        x = mat[2, 1] - mat[1, 2]
        y = mat[0, 2] - mat[2, 0]
        z = mat[1, 0] - mat[0, 1]
    vec = np.array([x, y, z])
    vec = vec / np.linalg.norm(vec)

    return vec * np.arccos(cosangle)


def rot_vec_to_axis_angle(rvec, angle_in_degree=True, allow_0_angle=True):
    """
    if allow_0_angle is False, then absolute 0 rotation vector will
    be treated as (0,0,2*pi); otherwise axis=(0,0,1) while angle=0
    """
    ang = np.linalg.norm(rvec)
    if ang > 0:
        angle = ang * 180.0 / np.pi if angle_in_degree else ang
        return (rvec / ang, angle)
    else:
        ang = ang if allow_0_angle else 360
        angle = ang if angle_in_degree else ang / 180.0 * np.pi
        return (np.array((0, 0, 1)), angle)


def vectors_to_qs(xyzs, keV, unit='', calc_angles=False):
    """
    from xyz lab coordinates of an outgoing beam vector, find q vectors based on beam energy.

    xyzs         :: should be an array with any shape, as long as the size of
                    the last dimension = 3
    keV          :: photon energy in keV
    unit = 'keV' :: output will be in unit of keV
         = 'test':: output will just be k_out_hat - k_in_hat, no units
           other :: output default in SI unit of m^-1
    calc_angles  :: if True, output will be a tuple of (q,tth,phi), where tth & phi
                    are the 2-theta scattering angle and the azimuthal angle relative
                    to the YZ plane;
                    if False, output will just be an array of q vectors
    """

    pxyzs = np.array(xyzs)

    ## get kout_hat
    pxyz_norms = np.sqrt(np.sum(pxyzs ** 2, axis=-1, keepdims=True))
    pxyzs = pxyzs / pxyz_norms

    k_in_hat = np.array([0.0, 0.0, 1.0])

    ## subtraction: k_out_hat - k_in_hat
    qs = pxyzs - k_in_hat

    ## now multiply by magnitude of k0 if needed
    if unit == 'keV':
        qs = qs * keV  # in unit of keV
    elif unit == 'test':
        pass  # check q directions only
    else:
        from scipy import constants  # in unit of m^-1
        qs = qs * (keV * 1.0e3 * constants.e / constants.hbar / constants.c)

    ## calc angles if needed
    if calc_angles:
        pxy_rhos = np.sqrt(np.sum(pxyzs[..., :2] ** 2, axis=-1))
        tths = np.arctan(pxy_rhos)
        phis = np.arcsin(abs(pxyzs[..., 0]) / np.fmax(pxy_rhos, 1.0e-8))  # 0 <= phi <= PI/2
        return qs, tths, phis
    else:
        return qs


class Detector(object):
    """
    Simple class that defines dimensions of an area detector
    """

    def __init__(self, npixel_x, npixel_y, pixelsize_x_mm, pixelsize_y_mm, desc=''):
        self.description = desc
        self.npx = npixel_x
        self.npy = npixel_y
        self.psizex = pixelsize_x_mm
        self.psizey = pixelsize_y_mm


class DetectorGeometry(object):
    """
    geometry of a detector as defined at 34IDE
    """

    def __init__(self, detname, det, trans_vec, rot_vec):
        self._rmatrix = None
        assert isinstance(det, Detector)
        self.detname = detname
        self.det = det
        self.translation = np.array(trans_vec).astype('float64')
        self.rotation = np.array(rot_vec).astype('float64')
        self.update_rmatrix()

    def update_rmatrix(self):
        """
        Call this after self.rotation is updated
        """
        self._rmatrix = rot_vec_to_matrix(self.rotation)

    @property
    def rmatrix(self):
        return np.copy(self._rmatrix)

    @property
    def inv_rmatrix(self):
        return self.rmatrix.T

    @property
    def proj_center(self):
        """
        return coordinates of the projected position of the sample on the
        detector plane
        """
        dx, dy = self.translation[:2]
        dpx = dx / self.det.psizex
        dpy = dy / self.det.psizey
        px = (self.det.npx - 1) * 0.5 - dpx
        py = (self.det.npy - 1) * 0.5 - dpy
        return px, py

    def pixels_to_xyzs(self, pixelIs, pixelJs):
        """
        return xyz coordiantes for each pixel, pixels' (i,j)-coordinates
        given by 2 separate arrays (or just 2 numbers)
        """
        pixIs = np.array(pixelIs)  # so the input can be either numbers or arrays of numbers
        pixJs = np.array(pixelJs)

        inputshape = pixIs.shape
        outputshape = tuple(list(inputshape) + [3])

        Zs = np.zeros(inputshape)
        pxyzs = np.column_stack((pixIs.flat, pixJs.flat, Zs.flat))
        pxyzs.shape = outputshape

        ## as if pixel (0,0) is at origin
        # default detector orientation: X outboard, Y up, in units of mm
        pxyzs[..., 0] *= self.det.psizex
        pxyzs[..., 1] *= self.det.psizey

        ## translate to center of detector, unit = mm
        pxyzs[..., 0] -= self.det.psizex * (self.det.npx - 1) * 0.5
        pxyzs[..., 1] -= self.det.psizey * (self.det.npy - 1) * 0.5

        ## translate acoording to geometry 
        pxyzs += self.translation

        ## now rotate according to geometry ##
        pxyzs = np.tensordot(pxyzs, self.rmatrix, axes=(-1, 1))

        return pxyzs

    def xyzs_to_pixels(self, xyzs):
        """
        for one/multiple given direction vector(s) xyzs, find
        the pixel indices that this direction would interset the
        detector plane.
        """
        inputshape = xyzs.shape
        assert inputshape[-1] == 3  # last dim should be 3

        nvectors = len(xyzs) // 3

        xyzs = np.reshape(xyzs, (nvectors, 3))  # convert to a Nx3 array for handling

        ### rotate xyzs back to non-rotated detector frame
        uxyzs = np.tensordot(xyzs, self.inv_rmatrix, axes=(-1, 1))

        ### find interesection pixel indices ###
        def get_pixel(uxyz):
            """
            To be called with map() to get pixel indices from input vector uxyz
            """
            ## find intersection of direction uxyz with plane z = z_dettrans
            ux, uy, uz = uxyz
            tx, ty, tz = self.translation

            if uz * tz <= 0:  # uxyz pointing away from detector plane
                px = py = np.nan
            else:  # px,py are the coordinates of intersection
                px = ux * tz / uz
                py = uy * tz / uz
            ## convert to pixel indices
            # x,y relative to center
            px -= tx
            py -= ty
            # x,y relative to (0,0)
            px += self.det.psizex * (self.det.npx - 1) * 0.5
            py += self.det.psizey * (self.det.npy - 1) * 0.5

            pxi = px / self.det.psizex
            pyi = py / self.det.psizey
            return np.array([pxi, pyi])

        pixels = np.array(list(map(get_pixel, uxyzs)))

        ### reshape output to oonform to the input shape, with last dim = 2
        outputshape = list(inputshape)
        outputshape[-1] = 2
        pixels.shape = tuple(outputshape)

        return pixels


# this function doesn't really belong to this module, but is no where else to 
#   put at this moment
def qs_to_ortho_hkls(qs, rec_latt_consts, samp_rot_matrix):
    """
    Convert q vectors to hkls in sample reciprocal space, for orthorgonal lattice systems only
    (cubic, tetragonal, orthorhombic)

    qs              :: should be an array with the size of the last dim = 3
    rec_latt_consts :: a 3-array storing the 3 reciprocal lattice constants a*,b*,c*
    * qs & rec_latt_consts should be of the same units
    samp_rot_matrix :: a 3x3 matrix that rotate sample orientation from standard to acutal

    returns q vectors in sample hkls (reciprocal lattice units)
    """
    rot_inverse = samp_rot_matrix.T  # the transverse of rotation matrix is its inverse

    ## mulitply the inverse rotation to the q vectors to get them in crystal coordinates
    qs_sample = qs.T  # reverses its dimensions from (n1,n2,...,3) to (3,...,n2,n1)
    qs_sample = np.tensordot(rot_inverse, qs_sample, axes=(1, 0))  # result of the dot is still of shape (3,...,n2,n1)
    qs_sample = qs_sample.T  # now reverse the dimensions back to (n1,n2,...,3)

    ## reduce the vectors by r.l.u.
    return qs_sample / rec_latt_consts
