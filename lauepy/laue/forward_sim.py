# -*- coding: utf-8 -*-
"""
Classes for simulating a Laue diffraction on 2D area detectors

classes:
    Reflection
    LaueSimulation

Author : Xu, Ruqing
Created: Thu Jan 14 18:36:04 2016
"""
import collections
from copy import deepcopy
from functools import wraps
from functools import partial as partialf
from time import perf_counter as timemeth

import numpy as np

from lauepy.rxlibs.xmd34 import geometry as geo
from lauepy.rxlibs.xmd34 import lattice as latt
from lauepy.rxlibs.xmd34 import utils


def measure(func):
    @wraps(func)
    def _time_it(*args, **kwargs):
        start = int(round(timemeth() * 1000))
        try:
            return func(*args, **kwargs)
        finally:
            end_ = int(round(timemeth() * 1000)) - start
            # print(f"Total execution time: {end_ if end_ > 0 else 0} ms")

    return _time_it


class Reflection(object):
    """
    Collection of related quanties & methods of a Bragg reflection

    Usage
    -----
    ::
        Reflection(hkl)

    Parameters
    ----------
        hkl : array-like
            an array with 3 elements that defines a reflection

    Methods
    -------
        gen_q_vec(xtal,rotation_matrix)
            generate q vector based on crystal lattice & rotation matrix
        calc_coords_on_detector(detector_geometry)
            compute pixel coordinates of this reflection according to the
            given detector_geometry

    Attibutes / Properties
    ----------------------
        hkl
            input hkl vector
        hkllabel
            (property) a string with the h,k,l numbers
        q
            scattering vector, calculated by calling self.gen_q_vec
        keV
            x-ray energy needed for this reflection, calculated by calling
            self.gen_q_vec
        keVlabel
            (property) a string with keV number printed
        k_out
            outgoing beam wavevector, calculated by calling self.gen_q_vec
        k_out_normalized
            (property) normalized k_out vector
        det_coords
            dict object generated by calling self.calc_coords_on_detector,
            has the form of {'detname1':(x1,y1),'detname2':(x2,y2),...}, with
            the keys being the "detname" filed of a given detector_geometry
        structure_factor
            defined for other programs to use, set to default of -1.0;
        sflabel
            (property) a string with the structure factor number
    """

    def __init__(self, hkl):
        self.hkl = np.array(hkl).astype('float32')
        self.q = None
        self.keV = None
        self.k_out = None
        self.det_coords = {}
        self.structure_factor = -1.0

    def gen_q_vec(self, xtal, rotation):
        """generate q vector based on xtal info and input rotation matrix"""
        self.q = np.dot(rotation, xtal.get_q_vec(self.hkl))
        self._calc_beam_k()

    def _calc_beam_k(self):
        """
        calculate beam energy and outgoing wavevector k_out
        vector, with option to be normalized
        """
        if self.q is not None:
            qx, qy, qz = self.q
            if qz < 0:  # only accessible if q is in the upstream half
                self.k_out = np.array([qx, qy, (-qz ** 2 + qx ** 2 + qy ** 2) / (-qz * 2)])
                k = (qx ** 2 + qy ** 2 + qz ** 2) / (-2 * qz)
                self.keV = utils.inv_m_to_keV(k)
            else:
                self.keV = np.nan

    @property
    def k_out_normalized(self):
        if self.k_out is not None:
            return self.k_out / np.linalg.norm(self.k_out)
        else:
            return None

    def calc_coords_on_detector(self, detgeometry):
        """
        Computes self.det_coords, which is a dict object with the form of
        {'detname1':(x1,y1),'detname2':(x2,y2),...}

        Input
        -----
            detgeometry
                should be of DetectorGeometry class
        """
        assert isinstance(detgeometry, geo.DetectorGeometry)
        detname = detgeometry.detname
        if self.k_out is None:
            self.det_coords[detname] = None
        else:
            coords = detgeometry.xyzs_to_pixels(self.k_out)
            if np.sum(np.isnan(coords)):
                self.det_coords[detname] = None
            else:
                self.det_coords[detname] = coords

    @property
    def hkllabel(self):
        return '{:g} {:g} {:g}'.format(self.hkl[0], self.hkl[1], self.hkl[2])

    @property
    def keVlabel(self):
        if self.keV:
            return '{:.3f}'.format(self.keV)
        else:
            return None

    @property
    def sflabel(self):
        return '{:.1f}'.format(self.structure_factor)


class LaueSimulation(object):
    """
    Collection of methods & data objects for a tweakable Laue Simulation

    Inputs:
        detectors & geometries (allows multiple)
        crystal structure
        crystal orientation matrix (optional)
        beam energy cutoff (optional)
        extra plot range (optional)
    Data objects:
        peakpool
        peaks to plot
    User Methods:
        sample rotate around X/Y/Z/H/F/particular pixel on detector
        detector rotate around X/Y/Z/beam-center
        detector translations
        get sample rotate matrix/vector
        get detector rotate vector
        close
    """

    def __init__(self, detgeos, xtal, xtalrotation=None, keV_cutoffs=(7, 30),
                 plot_range_padding_percent=0):
        self.plotslist = None
        self.hkls = None
        self.ys = None
        self.xs = None
        self.peaks_plot = None
        self.peakpool = None
        self.detgeos = []
        if isinstance(detgeos, geo.DetectorGeometry):
            self.detgeos.append(deepcopy(detgeos))
        elif isinstance(detgeos, collections.Sequence) and \
                isinstance(detgeos[0], geo.DetectorGeometry):
            self.detgeos.extend(detgeos)

        assert isinstance(xtal, latt.Xtal)
        self.xtal = xtal

        if xtalrotation is None:
            self.xtal_rot = np.identity(3)
        else:
            self.xtal_rot = np.array(xtalrotation)  # this will make an copy
            assert self.xtal_rot.shape == (3, 3)

        self.energy_cutoffs_keV = np.array(keV_cutoffs)

        self.plot_padding_ratio = plot_range_padding_percent / 100.0

        self.rotate_step_default = 5 * np.pi / 180.0  # default rotation step is 5 degrees

        # make copies of quanties that one might want to reset to after changes
        self._init_detgeos = deepcopy(self.detgeos)
        self._init_xtal_rot = deepcopy(self.xtal_rot)

        # default configurations for matplotlib
        self.plotsconfig = {'hkllabels': False,
                            'keVlabels': False,
                            'sflabels': False,
                            'fontsize': 20,
                            'marker': 'o',
                            'markersize': 81,
                            'color': 'r',
                            'facecolors': 'none',
                            'edgecolors': 'r'
                            }

    def reset_geometry(self):
        """
        Reset detector geometries and xtal rotation to initial values.
          re-generate peaks and plot if there are any
        """
        self.detgeos = deepcopy(self._init_detgeos)
        self.xtal_rot = deepcopy(self._init_xtal_rot)
        if self.peaks_plot is not None:
            self.gen_peaks_to_plot()
            self.update_plots()

    def det_index(self, detid):
        """
        Return index of a detector in from its name,
        if input is an integer, then check its range.
        Raises ValueError if input not valid.
        """
        ### check detector id ###
        wrongid = False
        if detid in range(len(self.detgeos)):
            idet = detid
        elif isinstance(detid, str):
            try:
                idet = [detgeo.detname for detgeo in self.detgeos].index(detid)
            except:
                wrongid = True
        else:
            wrongid = True

        if wrongid:
            raise ValueError('Unrecognized detector id: {}'.format(detid))

        return idet

    def gen_peaks_to_plot(self):
        """
        generate list of peaks to be plotted on each detector
        """
        if self.peakpool is None:
            #### generate a pool of peaks that can potentially hit detector ####

            energy_cutoffs_inv_m = utils.keV_to_inv_m(self.energy_cutoffs_keV)

            ### generate a bunch of hkls to be considered ###
            maxhkl = np.rint(energy_cutoffs_inv_m[1] / self.xtal.abc_stars * 2.1)
            hh = np.arange(-maxhkl[0], maxhkl[0] + 1)
            kk = np.arange(-maxhkl[1], maxhkl[1] + 1)
            ll = np.arange(-maxhkl[2], maxhkl[2] + 1)
            hs, ks, ls = np.meshgrid(hh, kk, ll, indexing='ij')
            hklpool = np.concatenate((hs[..., None], ks[..., None], ls[..., None]), axis=3)
            hklpool.shape = (len(hh) * len(kk) * len(ll), 3)

            ### check magnitude of q ###
            hklpool = [hkl for hkl in hklpool if
                       self.xtal.get_q_mag(hkl) <= energy_cutoffs_inv_m[1] * 2 and np.linalg.norm(hkl) > 0.1]

            ### check lattice symmetry selection ###
            hklpool = list(filter(self.xtal.reflection_is_allowed, hklpool))

            ### save these peaks as they can potentially appear  #
            #   on the detector through rotations              ###
            self.peakpool = [Reflection(hkl) for hkl in hklpool]

        #### if peak pool already exists, find those that will  #
        #    actually appear on each detector (or within the    #
        #    extended plotting ranges)                       ####

        ### further filter out hkls that's not included in the  #
        #   max Ewald sphere                                  ###
        peaks = list(filter(self.__peak_energy_verify, self.peakpool))

        ### get all peaks that are to hit the plot range ###
        npxys = np.array([(detgeo.det.npx, detgeo.det.npy) for
                          detgeo in self.detgeos])
        pr = self.plot_padding_ratio
        # plotranges[i] = ((xmin,xmax),(ymin,ymax)) for detector[i]
        plotranges = [((-pr * npxy[0], npxy[0] * (1 + pr)),
                       (-pr * npxy[1], npxy[1] * (1 + pr))) for npxy in npxys]

        detpeaks = []
        for i in range(len(self.detgeos)):
            # filter_f = lambda pk:peak_coords_calc_n_verify(i,pk)
            filter_f = partialf(self.__peak_coords_calc_n_verify, idet=i, plotrange=plotranges[i])
            detpeaks.append(list(filter(filter_f, peaks)))

        ### compute structure factors ###
        for i in range(len(self.detgeos)):
            for peak in detpeaks[i]:
                peak.structure_factor = self.xtal.get_structure_factor(peak.hkl)

        ### remove duplicates (higher harmonics) from list ###
        self.peaks_plot = []
        for i in range(len(self.detgeos)):
            filter_f = partialf(self.__not_a_harmonic, peaklist=detpeaks[i])
            # filter_f = lambda pk: self.__not_a_harmonic(pk,detpeaks[i])
            self.peaks_plot.append(list(filter(filter_f, detpeaks[i])))
        pks = self.peaks_plot[0]
        detname = self.detgeos[0].detname
        self.xs = [pk.det_coords[detname][0] for pk in pks]
        self.ys = [pk.det_coords[detname][1] for pk in pks]
        self.hkls = [pk.hkllabel for pk in pks]

    def _calc_structure_factors(self):
        pass

    def make_plots(self, axes_, **kwargs):
        """
        Make plot of simulated peaks on the given matplotlib Axes object
        (or a list of them for multiple detectors)
        axes_ must already have been generated by external program

        kwargs overides default configurations of plotting, check
        the dict of self.plotsconfig for available fields
        """
        for arg, val in list(kwargs.items()):
            self.plotsconfig[arg] = val

        if isinstance(axes_, collections.Sequence):
            n_axes = len(axes_)
            axeslist = list(axes_)
        else:
            n_axes = 1
            axeslist = [axes_]

        nplot = min(len(self.detgeos), n_axes)  # in case the numbers don't match

        self.plotslist = []

        for i in range(nplot):
            # do these for each plot:
            thisplot = {}
            detname = self.detgeos[i].detname

            ax = axeslist[i]
            thisplot['axes'] = ax

            pks = self.peaks_plot[i]
            self.xs = [pk.det_coords[detname][0] for pk in pks]
            self.ys = [pk.det_coords[detname][1] for pk in pks]

            ll = ax.scatter(self.xs, self.ys, marker=self.plotsconfig['marker'],
                            s=self.plotsconfig['markersize'],
                            facecolors=self.plotsconfig['facecolors'],
                            edgecolors=self.plotsconfig['edgecolors'])
            thisplot['points'] = ll

            labellist = []
            for pk in pks:
                if self.plotsconfig['hkllabels']:
                    ann = ax.annotate(pk.hkllabel, xy=pk.det_coords[detname],
                                      xytext=(+10, 6), textcoords='offset points',
                                      fontsize=self.plotsconfig['fontsize'],
                                      color=self.plotsconfig['color'])
                    labellist.append(ann)
                if self.plotsconfig['keVlabels']:
                    ann = ax.annotate(pk.keVlabel, xy=pk.det_coords[detname],
                                      xytext=(+10, -18), textcoords='offset points',
                                      fontsize=self.plotsconfig['fontsize'],
                                      color=self.plotsconfig['color'])
                    labellist.append(ann)
                if self.plotsconfig['sflabels']:
                    voffset = -36 if self.plotsconfig['keVlabels'] else -18
                    ann = ax.annotate(pk.sflabel, xy=pk.det_coords[detname],
                                      xytext=(-9, voffset), textcoords='offset points',
                                      fontsize=self.plotsconfig['fontsize'],
                                      color=self.plotsconfig['color'])
                    labellist.append(ann)
            thisplot['labels'] = labellist

            self.plotslist.append(thisplot)

    def update_plots(self, **kwargs):
        if self.plotslist is not None:
            axeslist = []
            for plot1 in self.plotslist:
                plot1['points'].remove()
                for label in plot1['labels']:
                    label.remove()
                axeslist.append(plot1['axes'])
            self.make_plots(axeslist, **kwargs)

    def clear_plots(self):
        if self.plotslist is not None:
            for plot1 in self.plotslist:
                plot1['points'].remove()
                for label in plot1['labels']:
                    label.remove()
        self.plotslist = None

    def rotate_xtal(self, axis, degree=None):
        """
        Method for user to rotate crystal and re-calcualte the simulation
        axis can be of different types:
            'X', 'Y', 'Z', 'H' or 'F': pre-defined beamline axis
            'pcenter' : projected center of the sample on the detector plane
            (x,y) : rotate around detector pixel given by (i,j)
            (u,v,w) : an arbitrary axis given by this vector
        degree : angle to be rotated in degree, optional
        """
        ### define axis ###
        if len(axis) == 2:
            # assume axis=(x,y) refers to a pixel on detector
            pixelxyz = self.detgeo.pixels_to_xyzs(*axis)
            rvec = geo.pxyzs_to_qs(pixelxyz, 1.0)
            rvec.shape = (3,)
        elif len(axis) == 3:
            rvec = np.array(axis)
        elif axis in 'Xx':
            rvec = np.array((1, 0, 0))
        elif axis in 'Yy':
            rvec = np.array((0, 1, 0))
        elif axis in 'Zz':
            rvec = np.array((0, 0, 1))
        elif axis in 'Hh':
            rvec = np.array((0, 1, 1))
        elif axis in 'Ff':
            rvec = np.array((0, -1, 1))
        elif axis == 'pcenter':
            px, py = self.detgeo.proj_center
            pixelxyz = self.detgeo.pixels_to_xyzs(px, py)
            rvec = geo.vectors_to_qs(pixelxyz, 1.0)
            rvec.shape = (3,)
        else:
            raise ValueError("Unrecognized input '{}' for 'axis'".format(axis))
        # normalize axis vector
        rvec = rvec / np.linalg.norm(rvec)

        ### multiply by rotation angle ###
        if degree is None:
            rvec *= self.rotate_step_default
        else:
            rvec *= degree / 180.0 * np.pi

        ### make new rotation matrix ###
        rmat = geo.rot_vec_to_matrix(rvec)
        self.xtal_rot = np.tensordot(rmat, self.xtal_rot, axes=(1, 0))

        ### update peak list ###
        if self.peaks_plot is not None:
            self.gen_peaks_to_plot()
            self.update_plots()

    def rotate_det(self, detid, axis, degree=None):
        """
        Method for user to rotate crystal and re-calcualte the simulation
        detid can either be an integer index or a string name of the detector.
        axis can be of different types:
            'X', 'Y', 'Z': pre-defined beamline axis
            'pcenter' : projected center of the sample on the detector plane
            (u,v,w) : an arbitrary axis given by this vector
        degree : angle to be rotated in degree, optional
        """
        idet = self.det_index(detid)
        detgeo = self.detgeos[idet]

        ### define axis ###
        if len(axis) == 3:
            rvec = np.array(axis)
        elif axis in 'Xx':
            rvec = np.array((1, 0, 0))
        elif axis in 'Yy':
            rvec = np.array((0, 1, 0))
        elif axis in 'Zz':
            rvec = np.array((0, 0, 1))
        elif axis == 'pcenter':
            # the proj. center rotation for the detector is different from that
            # for the xtal in that it's not around the q vector, but around the
            # position vector of the proj. center
            px, py = detgeo.proj_center
            rvec = detgeo.pixels_to_xyzs(px, py)
            rvec.shape = (3,)  # this line may be redundant
        else:
            raise ValueError("Unrecognized input '{}' for 'axis'".format(axis))
        # normalize axis vector
        rvec = rvec / np.linalg.norm(rvec)

        ### multiply by rotation angle ###
        if degree is None:
            rvec *= self.rotate_step_default
        else:
            rvec *= degree / 180.0 * np.pi
        ### make new rotation matrix ###
        rmat = geo.rot_vec_to_matrix(rvec)
        newdetrm = np.tensordot(rmat, detgeo.rmatrix, axes=(1, 0))
        newdetrv = geo.rot_matrix_to_vec(newdetrm)
        detgeo.rotation = newdetrv
        detgeo.update_rmatrix()

        ### update peak list ###
        if self.peaks_plot is not None:
            self.gen_peaks_to_plot()
            self.update_plots()

    def translate_det(self, detid, vec, distance=None):
        """
        User method to change detector translation vector and then update
        the calculated peaks.
        If distance is not given, then vec is used as absolute values;
        otherwise, the change will be made in the direction of vec with
        the specified distance.
        """
        idet = self.det_index(detid)

        if distance is None:
            distance = 1.0

        v = np.array(vec).astype('float32')
        v = v / np.linalg.norm(v) * distance

        self.detgeos[idet].translation += v

        ### update peak list ###
        if self.peaks_plot is not None:
            self.gen_peaks_to_plot()
            self.update_plots()

    ##########################################################
    #### internal util functions for class LaueSimulation ####
    ##########################################################
    def __peak_energy_verify(self, pk):
        """check if peak energy is smaller than max cutoff."""
        assert isinstance(pk, Reflection)
        pk.gen_q_vec(self.xtal, self.xtal_rot)
        if pk.keV < self.energy_cutoffs_keV[1]:
            return True
        else:
            return False

    @staticmethod
    def __not_a_harmonic(pk, peaklist):
        """
        check if a peak is a higher hamonic of any one in peaklist
        """
        ratio = np.ones(3)
        result = True
        ## loop over peaks, see if there's a shorter hkl  #
        #   that is also parallel to pk                  ##
        for peak in peaklist:
            parallel = None
            islonger = False
            for i in (0, 1, 2):
                a = pk.hkl[i]
                b = peak.hkl[i]
                if a == 0 and b == 0:
                    ratio[i] = np.nan
                elif b == 0 or a == 0:
                    # just one of a or b is zero, obviously not parallel
                    parallel = False
                    break
                elif abs(a) <= abs(b):
                    islonger = True
                    break
                else:
                    ratio[i] = a * 1.0 / b  # *1.0 only necessary in Py2

            if (parallel is False) or islonger:
                # if surely not parallel, or peak is longer in
                # any dimension, go to next peak
                continue
            ratiof = list(filter(np.isfinite, ratio))

            if len(ratiof) == 1:
                parallel = True
            elif len(ratiof) == 2:
                parallel = (ratiof[0] == ratiof[1])
            else:
                parallel = (ratio[0] == ratio[1]) and (ratio[0] == ratio[2])

            if parallel:
                result = False
                break
        return result

    def __peak_coords_calc_n_verify(self, pk, idet, plotrange):
        """
        Calculate pixel coordinates on i-th detector and
        check if the coordinates fall in the plot range
        """
        detgeo = self.detgeos[idet]

        pk.calc_coords_on_detector(detgeo)
        if pk.det_coords[detgeo.detname] is None:
            return False
        else:
            x, y = pk.det_coords[detgeo.detname]
            xmin, xmax = plotrange[0]
            ymin, ymax = plotrange[1]
            return x > xmin and x < xmax and y > ymin and y < ymax


class FastLaueSimulation_list(object):
    """
    Collection of methods & data objects for a tweakable Laue Simulation

    Inputs:
        detectors & geometries (allows multiple)
        crystal structure
        crystal orientation matrix (optional)
    Data objects:
        peakpool
    User Methods:
        detector rotate around X/Y/Z/beam-center
        detector translations
        get sample rotate matrix/vector
        get detector rotate vector
        close
    """

    def __init__(self, hkl, detgeos, xtal, xtalrotation=None):
        self.detgeos = []
        if isinstance(detgeos, geo.DetectorGeometry):
            self.detgeos.append(deepcopy(detgeos))
        elif isinstance(detgeos, collections.Sequence) and \
                isinstance(detgeos[0], geo.DetectorGeometry):
            self.detgeos.extend(detgeos)

        # self.hkl=eulerparsehkl.extract_hkl(filename)
        self.hkl = hkl
        self.peakpool = [Reflection(hkl) for hkl in self.hkl]

        assert isinstance(xtal, latt.Xtal)
        self.xtal = xtal

        if xtalrotation is None:
            self.xtal_rot = np.identity(3)
        else:
            self.xtal_rot = np.array(xtalrotation)  # this will make an copy
            assert self.xtal_rot.shape == (3, 3)

        self.rotate_step_default = 5 * np.pi / 180.0  # default rotation step is 5 degrees

        # make copies of quanties that one might want to reset to after changes
        self._init_detgeos = deepcopy(self.detgeos)
        self._init_xtal_rot = deepcopy(self.xtal_rot)

    def reset_geometry(self):
        """
        Reset detector geometries and xtal rotation to initial values.
          re-generate peaks and plot if there are any
        """
        self.detgeos = deepcopy(self._init_detgeos)
        self.xtal_rot = deepcopy(self._init_xtal_rot)
        self.gen_det_peaks()

    def det_index(self, detid):
        """
        Return index of a detector in from its name,
        if input is an integer, then check its range.
        Raises ValueError if input not valid.
        """
        ### check detector id ###
        wrongid = False
        if detid in range(len(self.detgeos)):
            idet = detid
        elif isinstance(detid, str):
            try:
                idet = [detgeo.detname for detgeo in self.detgeos].index(detid)
            except:
                wrongid = True
        else:
            wrongid = True

        if wrongid:
            raise ValueError('Unrecognized detector id: {}'.format(detid))

        return idet

    @measure
    def gen_det_peaks(self):

        # The old self.__peak_coords_calc_n_verify modified the peaks to actually compute the position on the detector.
        # So need to lift that from the old method assuming a single detector.  Might delete multidet capability.
        detgeo = self.detgeos[0]
        detname = detgeo.detname
        for pk in self.peakpool:
            pk.gen_q_vec(self.xtal, self.xtal_rot)
            pk.calc_coords_on_detector(detgeo)
        #             print(np.round(pk.keV,0),'keV')
        self.xs = [pk.det_coords[detname][0] for pk in self.peakpool]
        self.ys = [pk.det_coords[detname][1] for pk in self.peakpool]
        self.hkls = [pk.hkllabel for pk in self.peakpool]

    def rotate_det(self, detid, axis, degree=None):
        """
        Method for user to rotate crystal and re-calcualte the simulation
        detid can either be an integer index or a string name of the detector.
        axis can be of different types:
            'X', 'Y', 'Z': pre-defined beamline axis
            'pcenter' : projected center of the sample on the detector plane
            (u,v,w) : an arbitrary axis given by this vector
        degree : angle to be rotated in degree, optional
        """
        idet = self.det_index(detid)
        detgeo = self.detgeos[idet]

        ### define axis ###
        if len(axis) == 3:
            rvec = np.array(axis)
        elif axis in 'Xx':
            rvec = np.array((1, 0, 0))
        elif axis in 'Yy':
            rvec = np.array((0, 1, 0))
        elif axis in 'Zz':
            rvec = np.array((0, 0, 1))
        elif axis == 'pcenter':
            # the proj. center rotation for the detector is different from that
            # for the xtal in that it's not around the q vector, but around the
            # position vector of the proj. center
            px, py = detgeo.proj_center
            rvec = detgeo.pixels_to_xyzs(px, py)
            rvec.shape = (3,)  # this line may be redundant
        else:
            raise ValueError("Unrecognized input '{}' for 'axis'".format(axis))
        # normalize axis vector
        rvec = rvec / np.linalg.norm(rvec)

        ### multiply by rotation angle ###
        if degree is None:
            rvec *= self.rotate_step_default
        else:
            rvec *= degree / 180.0 * np.pi
        ### make new rotation matrix ###
        rmat = geo.rot_vec_to_matrix(rvec)
        newdetrm = np.tensordot(rmat, detgeo.rmatrix, axes=(1, 0))
        newdetrv = geo.rot_matrix_to_vec(newdetrm)
        detgeo.rotation = newdetrv
        detgeo.update_rmatrix()

        ### update peak list ###
        self.gen_det_peaks()

    def translate_det(self, detid, vec, distance=None):
        """
        User method to change detector translation vector and then update
        the calculated peaks.
        If distance is not given, then vec is used as absolute values;
        otherwise, the change will be made in the direction of vec with
        the specified distance.
        """
        idet = self.det_index(detid)

        if distance is None:
            distance = 1.0

        v = np.array(vec).astype('float32')
        v = v / np.linalg.norm(v) * distance

        self.detgeos[idet].translation += v
        self.gen_det_peaks()
