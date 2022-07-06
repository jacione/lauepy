"""
Written by Matt Wilkin, Yueheng Zhang, and Nick Porter
"""
# Laue Simulation and Deconvolution

import copy
import json
import random
import re
import subprocess as sub
import sys

import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from matplotlib import pyplot as plt
from progressbar import progressbar as pbar

import lauepy.laue.forward_sim as fsim
import lauepy.laue.peaks as pk
import lauepy.laue.rlv_to_spec as so
from lauepy.laue.disorientation import calc_disorient, rmat_2_quat
from lauepy.laue.pflibs import extract_rlv, rot_wcha_modified
from lauepy.rxlibs.xmd34 import geometry as geo
from lauepy.rxlibs.xmd34 import lattice as latt

find_eul = re.compile(r"EulerAngles\d\s+{\s*(.+)}")
find_gs = re.compile(r"\[\s+\d+]\s+\((.+?)\)")
find_rms = re.compile(r"rms_error\d\s+(.+?)\t")
find_goodness = re.compile(r"[$]goodness\d\t+(\d+.*\d+)\t")
find_cols = re.compile(r"{(.+?)}")
find_rlv = re.compile(r"recip_lattice\d\s+{(.+)}")
find_rmat = re.compile(r"[$]rotation_matrix\d\t+{(.+)}")
find_hkls = re.compile(r"\[\s+\d+]\s+\(.+?\)\s+\((.+?)\)")


class AutoLaue:

    def __init__(self, config):
        # These variables must be defined in self.index()
        self.config = config
        self.system = sys.platform
        self.working_dir = config['working_dir']
        self.peak_dict = pk.load_peaks(config)
        self.phichitheta = self.peak_dict['info']['angles']

        self.times = None
        self.comb_sub = None
        self.goodness = None
        self.pattern_ID = None
        self.pattern_dict = None
        self.tolerance = None
        self.frequency = None
        self.mis_err = None
        self.material = None
        self.space_group = None
        self.lattice_params = None
        self.atom_list = None
        self.xtal = None

        self.pix_x, self.pix_y = config['det_pixels']
        self.pitch_x, self.pitch_y = config['det_pitch']
        self.det_rotation = np.array(config['det_rotation'])
        self.det_translation = np.array(config['det_translation'])
        self.ccd1 = geo.Detector(self.pix_x, self.pix_y, self.pitch_x, self.pitch_y, config['det_name'])
        self.geo_ccd = geo.DetectorGeometry('ccd1', self.ccd1, self.det_translation, self.det_rotation)

    def set_params(self, substrate=False):
        name = 'sample'
        if substrate:
            name = 'substrate'

        self.times = self.config[f'laue_{name}_times']
        self.comb_sub = self.config[f'laue_{name}_comb_sub']
        self.goodness = self.config[f'laue_{name}_goodness']
        self.tolerance = self.config[f'laue_{name}_tolerance']
        self.frequency = self.config[f'laue_{name}_frequency']
        self.mis_err = self.config[f'laue_{name}_mis_err']
        with open(f'{self.config["lauepy_dir"]}/crystals/{self.config[name]}.json') as f:
            xtal_params = json.load(f)

        self.material = xtal_params['material']
        self.space_group = xtal_params['space_group']
        a, b, c, ang1, ang2, ang3 = xtal_params['lattice_params']
        self.lattice_params = f"{{ {a*1e9}, {b*1e9}, {c*1e9}, {ang1}, {ang2}, {ang3} }}"
        self.atom_list = [latt.AtomInCell(atom_pos[0], *atom_pos[1]) for atom_pos in xtal_params['pos_list']]
        self.xtal = latt.Xtal(a, b, c, ang1, ang2, ang3, atomlist=self.atom_list)

    def raster_beam(self, microstructure, stepx, stepy, beam_width, beam_height, num_peaks, xs, ys):
        w = beam_width
        h = beam_height
        beam_rast = []
        stack = []
        locations = []
        grain_list = []
        count = 0
        for i in range(0, microstructure.shape[0] - (h - stepy), stepy):
            for j in range(0, microstructure.shape[1] - (w - stepx), stepx):
                mic = copy.copy(microstructure)
                mic[i:i + h, j:j + w] *= 15
                beam_rast.append(mic)
                box = copy.copy(microstructure[i:i + h, j:j + w])
                grains = np.unique(box.ravel())
                grain_list.append([grain for grain in grains])

                img = np.zeros((self.pix_x, self.pix_y)).astype(np.uint16)
                locations.append([i, j])
                for l in range(len(grains)):

                    new_xs = xs[int(grains[l]) - 1, :num_peaks]
                    new_ys = ys[int(grains[l]) - 1, :num_peaks]

                    for r in range(num_peaks):
                        img[int(new_xs[r]), int(new_ys[r])] = 1e4
                img = gf(img, 2)

                stack.append(img)
                #                 print(count)
                count += 1
        stack = np.moveaxis(np.array(stack), 0, -1)

        return stack, beam_rast, grain_list, locations

    def calc_gs(self):
        # TODO This function would be a lot more transparent if it were vectorized.
        theta = np.sqrt(self.det_rotation[0]**2 + self.det_rotation[1]**2 + self.det_rotation[2]**2)
        # TODO theta = np.linalg.norm(self.det_rotation)
        c = np.cos(theta)
        s = np.sin(theta)
        c1 = 1 - c
        Rx, Ry, Rz = self.det_rotation / theta
        rho00 = c + Rx * Rx * c1
        rho01 = Rx * Ry * c1 - Rz * s
        rho02 = Ry * s + Rx * Rz * c1
        rho10 = Rz * s + Rx * Ry * c1
        rho11 = c + Ry * Ry * c1
        rho12 = -Rx * s + Ry * Rz * c1
        rho20 = -Ry * s + Rx * Rz * c1
        rho21 = Rx * s + Ry * Rz * c1
        rho22 = c + Rz * Rz * c1
        for frame_data in self.peak_dict.values():
            if 'coords' not in frame_data.keys():
                continue
            frame_data['G_vectors'] = [None for _ in frame_data['coords']]
            for i, coords in enumerate(frame_data['coords']):
                px, py = coords
                # start calculating the g-vectors below ###################################
                xd = (px - 0.5 * (self.pix_x - 1)) * self.pitch_x
                yd = (py - 0.5 * (self.pix_y - 1)) * self.pitch_y
                zd = 0
                # translate (xd,yd) by the vector P to get (xd,yd,zd)
                x, y, z = xd + self.det_translation[0], yd + self.det_translation[1], zd + self.det_translation[2]
                # rotate about R
                X = rho00 * x + rho01 * y + rho02 * z
                Y = rho10 * x + rho11 * y + rho12 * z
                Z = rho20 * x + rho21 * y + rho22 * z
                # normalize
                total = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
                # TODO total = np.linalg.norm([X,Y,Z])
                x_lab, y_lab, z_lab = X / total, Y / total, Z / total

                # normalize qhat
                q_x, q_y, q_z = x_lab, y_lab, z_lab - 1
                total_qhat = np.sqrt(q_x ** 2 + q_y ** 2 + q_z ** 2)
                # TODO total_qhat = np.linalg.norm([q_x, q_y, q_z])
                q_x = q_x / total_qhat
                q_y = q_y / total_qhat
                q_z = q_z / total_qhat

                frame_data['G_vectors'][i] = (q_x, q_y, q_z)
        pk.save_peaks(self.config, self.peak_dict)
        return

    def write_mat_params(self, gvectors):
        s = f"$filetype	PeaksFile\n" \
            f"$structureDesc {self.material}\n" \
            f"$latticeParameters {self.lattice_params}\n" \
            f"$latticeAlphaT {self.lattice_params[0]*10}E-10\n" \
            f"$lengthUnit nm\n" \
            f"$SpaceGroup {self.space_group}\n" \
            f"$N_Ghat+Intens {len(gvectors)}\n"
        for G in gvectors:
            s += f"{G[0]},{G[1]},{G[2]},1.0\n"
        with open(f'{self.working_dir}/peaks/Peaks.txt', 'w') as f:
            f.write(s)
        return

    def run_euler(self):
        # TODO: put the euler program for linux/mac/windows in the same directory and select based on current OS

        sub.run([
            f'{self.config["lauepy_dir"]}/scripts/eulerlinux',
            '-k', '16',
            '-t', '24',
            '-c', '72',
            '-a', f'{self.tolerance}',
            '-f', f'{self.working_dir}/peaks/Peaks.txt',
            '-o', f'{self.working_dir}/peaks/Index.txt',
            ],
            stdout=sub.DEVNULL,
            stderr=sub.DEVNULL
        )
        return

    @staticmethod
    def indx(g, known_gs):
        x = [np.degrees(np.arccos(np.array(g) @ np.array(val))) for val in known_gs]
        ind = np.argmin(x)

        return ind

    def extract_vals(self):
        with open(f"{self.working_dir}/peaks/Index.txt") as f:
            file_split = re.split("[$]pattern\d", f.read())[1:]
        found_patts = []

        for patt in file_split[:6]:

            goodness = np.round(float(find_goodness.findall(patt)[0]), 2)

            if goodness > self.goodness:
                cols = find_cols.findall(find_rlv.findall(patt)[0])
                vals = [m.split(',') for m in cols]
                rlv = np.array([[float(v) for v in x] for x in vals]).T
                rms = float(find_rms.findall(patt)[0])

                hkls = find_hkls.findall(patt)
                hkls = [[int(k) for k in hkl.split()] for hkl in hkls]

                cols = find_cols.findall(find_rmat.findall(patt)[0])
                vals = [m.split(',') for m in cols]
                rmat = np.array([[float(v) for v in x] for x in vals])

                gs = find_gs.findall(patt)
                gs = [[float(k) for k in g.split()] for g in gs]

                found_patts.append([hkls, gs, rmat, rlv, goodness, rms])

        return found_patts

    def orientation_counter(self, frame_data):
        peak_coords = np.array(frame_data['coords'])
        num_peaks = peak_coords.shape[0]

        # FIXME are these just magic numbers, or are they physically determined?
        if (num_peaks - 10) <= self.comb_sub:
            r = num_peaks
            size = 1
        elif num_peaks < 7:
            r = num_peaks
            size = 1
        else:
            r = num_peaks - self.comb_sub
            size = self.times

        indices = [random.sample(range(num_peaks), r) for _ in range(size)]

        reduced_list = []
        for count_index, ind in enumerate(indices):
            gvectors = np.array(frame_data['G_vectors'])[ind]

            self.write_mat_params(gvectors)
            self.run_euler()

            found_patts = self.extract_vals()
            ### 
            for patts in found_patts:
                hkls, gs, rmat, rlv, goodness, rms = patts
                phi, chi, theta = self.phichitheta
                spec_ori = so.transform(rlv.T, phi, chi, theta)
                ids = np.array([self.indx(g, gvectors) for g in gs])
                rtm1 = rmat

                dist_loss, xys = self.loss_function_distance(rmat.T, peak_coords[ids], hkls)

                No = True
                if len(reduced_list) != 0:
                    q1s = np.array([rtm1 for _ in reduced_list])

                    q2s = np.array([r[4] for r in reduced_list])

                    misorientations = calc_disorient(rmat_2_quat(q1s), rmat_2_quat(q2s))
                    for rr, r in enumerate(reduced_list):
                        rtm2 = r[4]

                        mis = misorientations[rr]

                        if mis < self.mis_err:
                            r[-1] += 1
                            r[0].append(rms)
                            r[1].append(goodness)
                            r[2].append(dist_loss)
                            r[3].append(len(gs))
                            r[-2].append(xys)

                            No = False
                            break
                if No:
                    reduced_list.append([[rms], [goodness], [dist_loss], [len(gs)], rmat, spec_ori, [xys], 1])

        confidence = np.array([r[-1] for r in reduced_list])
        good_ind = np.where(confidence >= self.frequency)
        reduced_list = [reduced_list[i] for i in list(good_ind[0])]
        return reduced_list

    def index(self):
        print("\nIndexing Laue peaks...")
        self.pattern_ID = 1
        self.pattern_dict = {}

        self.calc_gs()

        for frame, frame_data in pbar(self.peak_dict.items()):
            if frame.startswith('frame_'):
                i = int(frame[-5:])  # The frame number should be the last 5 characters in the dictionary key
                if frame_data['num_peaks'] < 3:
                    continue
            elif frame == 'substrate':
                i = -1
                self.set_params(substrate=True)
            else:
                continue
            reduced_o_list = self.orientation_counter(frame_data)

            for patt in reduced_o_list:
                rms, goodness, dist_loss, num_peaks, rmat, spec_ori, xys, conf = patt
                blah = np.array([len(xys) for _ in xys])

                self.pattern_dict[f'pattern_{self.pattern_ID}'] = {
                    'Rot_mat': rmat.tolist(),
                    'Goodness': float(np.mean(goodness)),
                    'Dist': float(np.mean(dist_loss)),
                    'Spec_Orientation': [s.tolist() for s in spec_ori],
                    'Center_Frame': i,
                    'Count': conf,
                    'RMS': float(np.mean(rms)),
                    'Num_Peaks': float(np.mean(num_peaks)),
                    'Pos': frame_data['positions'],
                    'xys': xys[np.argmax(blah)]
                }

                self.pattern_ID += 1
            if frame == 'substrate':
                self.set_params(substrate=False)
        with open(f"{self.working_dir}/peaks/patterns.json", 'w') as f:
            json.dump(self.pattern_dict, f)
        pk.save_peaks(self.config, self.peak_dict)
        return

    def loss_function_distance(self, xtal_rmat, xy_exp, hkl_labels):

        # x is detector parameter [P0, P1, P2, R0, R1, R2]
        # return the average distance between the measurement and the simulated the data
        #         tvec = np.array([x[0], x[1], x[2]])
        #         rvec = np.array([x[3], x[4], x[5]])
        # self.geoccd = geo.DetectorGeometry('ccd1', ccd1, tvec, rvec)
        #         hkl_labels = self.extract_hkl("Index.txt")
        # print("hkl_labels",hkl_labels)
        # print("self.geo_ccd",self.geo_ccd)
        # print("self.xtal",self.xtal)
        # print("xtal_rmat",xtal_rmat)
        xtalsimu = fsim.FastLaueSimulation_list(hkl_labels, self.geo_ccd, self.xtal, xtal_rmat)
        xtalsimu.gen_det_peaks()
        peak_list_r = np.array(xy_exp)
        # print(peak_list_r)
        peak_list_f = [xy for xy in zip(xtalsimu.xs, xtalsimu.ys)]

        d = cdist(peak_list_f, peak_list_r, metric='euclidean')
        indx_min = np.argmin(d, axis=1)
        distance_sum = 0
        for k4 in range(len(d)):
            distance_sum += d[k4][indx_min[k4]]

        # print("current R and P is", rvec, tvec)
        # print("average distance is", distance_sum / len(d))

        return distance_sum / len(d), peak_list_f

    @staticmethod
    def laue_transform(phi, chi, theta):
        data = extract_rlv()
        Mr = rot_wcha_modified(data, ['y', 'z', 'x'], [-theta, chi - 90, -phi])
        in_plane = Mr @ np.array([1, 0, 0])
        out_plane = Mr @ np.array([0, 1, 0])

        return in_plane, out_plane

    @staticmethod
    def rot_wcha_modified(data, axis, angle):
        for n in range(len(axis)):
            r = Rotation.from_euler(axis[n], angle[n], degrees=True).as_matrix()
            data = data @ r.T
        return data

    @staticmethod
    def extract_rlv():
        file_name = open("./Index.txt").read()
        file_split = re.split("[$]pattern\d", file_name)[1:]

        rlv_list = []
        for patt in file_split:
            rlv_string = find_rlv.findall(patt)[0]
            rlv_colum_str = re.findall(r"{[-+]?\d+\.?\d*,[-+]?\d+\.?\d*,[-+]?\d+\.?\d*}", rlv_string)

            rlv_colum_list = []
            for rlv_colum in rlv_colum_str:
                rlv_colum_list.append([float(s) for s in re.findall(r"[-+]?\d+\.?\d*", rlv_colum)])
            rlv_colum_array = np.array(rlv_colum_list)
            rlv_list.append(rlv_colum_array)

        return np.array(rlv_list[0])
