"""
Written by Matt Wilkin and Yueheng Zhang
"""
# Laue Simulation and Deconvolution

import copy
import json
import random
import re
import subprocess as sub

import numpy as np
from lauepy.rxlibs.xmd34 import geometry as geo
from lauepy.rxlibs.xmd34 import lattice as latt
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

import lauepy.laue.rlv_to_spec as so
import lauepy.laue.forward_sim as fsim
from lauepy.laue.disorientation import calc_disorient, rmat_2_quat
from lauepy.laue.pflibs import extract_rlv, rot_wcha_modified


class AutoLaue:

    def __init__(self, crystal_path, det_path, sym='Cubic_Sym.npy'):
        # These variables must be defined in self.index()
        self.params = None
        self.times = None
        self.comb_sub = None
        self.goodness = None
        self.system = None
        self.pattern_ID = None
        self.pattern_dict = None
        self.peak_dict = None
        self.group_dict = None
        self.tolerance = None
        self.frequency = None
        self.mis_err = None

        with open(crystal_path) as f:
            crystal_params = json.load(f)
        with open(det_path) as f:
            det_params = json.load(f)
        self.sym = sym

        self.find_eul = re.compile(r"EulerAngles\d\s+{\s*(.+)}")
        self.find_gs = re.compile(r"\[\s+\d+]\s+\((.+?)\)")
        self.find_rms = re.compile(r"rms_error\d\s+(.+?)\t")
        self.find_goodness = re.compile(r"[$]goodness\d\t+(\d+.*\d+)\t")
        self.find_cols = re.compile(r"{(.+?)}")
        self.find_rlv = re.compile(r"recip_lattice\d\s+{(.+)}")
        self.find_cols = re.compile(r"{(.+?)}")
        self.find_rmat = re.compile(r"[$]rotation_matrix\d\t+{(.+)}")
        self.find_hkls = re.compile(r"\[\s+\d+]\s+\(.+?\)\s+\((.+?)\)")
        pos_list = crystal_params['pos_list']
        self.material = crystal_params['material']

        self.space_group = crystal_params['space_group']
        a, b, c, ang1, ang2, ang3 = crystal_params['lattice_params']
        self.latticeParameters = "{ %s, %s, %s, %s, %s, %s }" % (a * 1e9, b * 1e9, c * 1e9, ang1, ang2, ang3)
        # self.atom_list = [latt.AtomInCell(self.material,*atom_pos) for atom_pos in pos_list]
        # atom_list = [latt.AtomInCell('Al',*atom_pos) for atom_pos in Al_pos_list] +\
        #             [latt.AtomInCell('O',*atom_pos) for atom_pos in O_pos_list]
        self.atom_list = [latt.AtomInCell(atom_pos[0], *atom_pos[1]) for atom_pos in pos_list]
        self.xtal = latt.Xtal(a, b, c, ang1, ang2, ang3, atomlist=self.atom_list)

        self.trans_vec = det_params['trans_vec']
        self.pix_x, self.pix_y, self.pitch_x, self.pitch_y, name = det_params['pixels']

        self.rot_vec = np.array(det_params['rot_vec'])
        self.trans_vec = np.array(det_params['trans_vec'])
        args = det_params['pixels']
        self.ccd1 = geo.Detector(*args)
        self.geo_ccd = geo.DetectorGeometry('ccd1', self.ccd1, self.trans_vec, self.rot_vec)
        self.det_params = det_params['pixels']

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
        theta = np.sqrt(
            self.rot_vec[0] * self.rot_vec[0] + self.rot_vec[1] * self.rot_vec[1] + self.rot_vec[2] * self.rot_vec[2])
        c = np.cos(theta)
        s = np.sin(theta)
        c1 = 1 - c
        Rx = self.rot_vec[0]
        Ry = self.rot_vec[1]
        Rz = self.rot_vec[2]
        Rx /= theta
        Ry /= theta
        Rz /= theta
        # print(Rx,Ry,Rz)
        rho00 = c + Rx * Rx * c1
        rho01 = Rx * Ry * c1 - Rz * s
        rho02 = Ry * s + Rx * Rz * c1
        rho10 = Rz * s + Rx * Ry * c1
        rho11 = c + Ry * Ry * c1
        rho12 = -Rx * s + Ry * Rz * c1
        rho20 = -Ry * s + Rx * Rz * c1
        rho21 = Rx * s + Ry * Rz * c1
        rho22 = c + Rz * Rz * c1
        for peak in self.peak_dict:
            peak_xy = self.peak_dict[peak]['XY']
            # start calculating the g-vectors below ###################################
            px = peak_xy[0]
            py = peak_xy[1]

            #             sizeX = 28.2/516
            #             sizeY = 28.2/517
            #             xd = (px - 0.5*(self.pix_x -1))*sizeX
            #             yd = (py - 0.5*(self.pix_y -1))*sizeY
            xd = (px - 0.5 * (self.pix_x - 1)) * self.pitch_x
            yd = (py - 0.5 * (self.pix_y - 1)) * self.pitch_y
            zd = 0
            # translate (xd,yd) by the vector P to get (xd,yd,zd)
            x, y, z = xd + self.trans_vec[0], yd + self.trans_vec[1], zd + self.trans_vec[2]
            # rotate about R
            X = rho00 * x + rho01 * y + rho02 * z
            Y = rho10 * x + rho11 * y + rho12 * z
            Z = rho20 * x + rho21 * y + rho22 * z
            # normalize
            total = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
            x_lab, y_lab, z_lab = X / total, Y / total, Z / total

            # normalize qhat
            q_x, q_y, q_z = x_lab, y_lab, z_lab - 1
            total_qhat = np.sqrt(q_x ** 2 + q_y ** 2 + q_z ** 2)
            q_x = q_x / total_qhat
            q_y = q_y / total_qhat
            q_z = q_z / total_qhat

            self.peak_dict[peak]['G_Vector'] = (q_x, q_y, q_z)

        return

    def write_mat_params(self, intensities, gvectors, ID):
        peak_file = open('Peaks.txt', 'w')
        lines = open("trial_gold_Peak.txt").readlines()
        i = 1
        while i <= 12:
            for s in lines:
                if s.startswith("$filetype") or \
                        s.startswith("// parameters defining the crystal structure") or \
                        s.startswith("$latticeAlphaT") or \
                        s.startswith("$lengthUnit") or \
                        s.startswith("// the following table contains xyz compotnents of G^ and the integral of the "
                                     "peak"):
                    peak_file.write(s)
                if s.startswith("$structureDesc"):
                    peak_file.write("$structureDesc		%s\n" % self.material)
                if s.startswith("$latticeParameters"):
                    peak_file.write("$latticeParameters	%s	// 2006, CODATA\n" % self.latticeParameters)
                if s.startswith("$SpaceGroup"):
                    peak_file.write(
                        "$SpaceGroup			%s					// Space Group number from International\n" % self.space_group)
                if s.startswith("$N_Ghat+Intens"):
                    peak_file.write("$N_Ghat+Intens 	%s		// number of G^ vectors\n" % len(intensities))
                i += 1
            # print("gvectors_group is",gvectors_group)
            for n, new_Value in enumerate(gvectors):
                peak_file.write("%s,%s,%s,%s\n" % (new_Value[0], new_Value[1], new_Value[2], intensities[n]))
            peak_file.close()

        return

    def run_euler(self):
        if self.system == 'linux':
            sub.call(
                ['./eulerlinux', '-k', '24', '-t', '24', '-c', '72', '-a', '%s' % self.tolerance, '-f', 'Peaks.txt'])
        elif self.system == 'mac':
            sub.call(['./euler', '-k', '30', '-t', '30', '-c', '180', '-f', 'Peaks.txt'])  ## for sapphire
            # sub.call(['./Euler','-k','24','-t','24','-c','72','-a','%s'%self.tolerance,'-f','Peaks.txt'])
        elif self.system == 'windows':
            sub.call(['euler', '-k', '16', '-t', '24', '-c', '72', '-a', '%s' % self.tolerance, '-f', 'Peaks.txt'])
        return

    @staticmethod
    def indx(g, known_gs):
        x = [np.degrees(np.arccos(np.array(g) @ np.array(val))) for val in known_gs]
        ind = np.argmin(x)

        return ind

    def extract_vals(self):
        file_name = open("Index.txt").read()
        file_split = re.split("[$]pattern\d", file_name)[1:]
        found_patts = []

        for patt in file_split[:6]:

            goodness = np.round(float(self.find_goodness.findall(patt)[0]), 2)

            if goodness > self.goodness:
                cols = self.find_cols.findall(self.find_rlv.findall(patt)[0])
                vals = [m.split(',') for m in cols]
                rlv = np.array([[float(v) for v in x] for x in vals]).T
                rms = float(self.find_rms.findall(patt)[0])

                hkls = self.find_hkls.findall(patt)
                hkls = [[int(k) for k in hkl.split()] for hkl in hkls]

                cols = self.find_cols.findall(self.find_rmat.findall(patt)[0])
                vals = [m.split(',') for m in cols]
                rmat = np.array([[float(v) for v in x] for x in vals])

                gs = self.find_gs.findall(patt)
                gs = [[float(k) for k in g.split()] for g in gs]

                found_patts.append([hkls, gs, rmat, rlv, goodness, rms])

        return found_patts

    def orientation_counter(self, group):
        all_peakIDs = self.group_dict[group]['ID_List']
        n = len(all_peakIDs)

        if (n - 10) <= self.comb_sub:
            r = n
            size = 1
        elif n < 7:
            r = n
            size = 1
        else:
            r = n - self.comb_sub
            size = self.times

        indices = [random.sample(range(n), r) for time in range(size)]
        print("%s combinations containing %s peaks" % (size, r))

        reduced_list = []
        count_index = 0
        for ind in indices:
            peakIDs = [all_peakIDs[i] for i in ind]

            gvectors = [self.peak_dict['peak_%s' % i]['G_Vector'] for i in peakIDs]

            self.write_mat_params(np.ones(len(gvectors)), gvectors, count_index)
            count_index += 1
            self.run_euler()

            found_patts = self.extract_vals()
            ### 
            for patts in found_patts:
                hkls, gs, rmat, rlv, goodness, rms = patts
                phi, chi, theta = self.params
                spec_ori = so.transform(rlv, phi, chi, theta)
                ids = np.array([peakIDs[self.indx(g, gvectors)] for g in gs])
                rtm1 = rmat

                dist_loss, xys = self.loss_function_distance(rmat.T,
                                                             [self.peak_dict['peak_%s' % ID]['XY'] for ID in ids], hkls)

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

    def index(self, param_path, substrate=False):
        with open(param_path) as f:
            params = json.load(f)
        self.times = params['times']
        self.comb_sub = params['comb_sub']
        self.goodness = params['goodness']
        self.system = params['system']

        self.pattern_ID = 1
        self.pattern_dict = {}
        with open(params['peak_dict']) as f:
            self.peak_dict = json.load(f)
        with open(params['group_dict']) as f:
            self.group_dict = json.load(f)

        self.tolerance = params['tolerance']
        self.frequency = params['frequency']
        self.mis_err = params['mis_err']
        self.calc_gs()
        count = 0

        for group in self.group_dict:
            print('Frame', self.group_dict[group]['Center_Frame'])
            count += 1
            self.params = self.group_dict[group]['phichitheta']
            reduced_o_list = self.orientation_counter(group)

            for patt in reduced_o_list:
                rms, goodness, dist_loss, num_peaks, rmat, spec_ori, xys, conf = patt
                #                 ids = list(ids)
                #                 for i in ids:
                #                     self.peak_dict['peak_%s'%i]['Pattern_ID'] = 'pattern_%s'%self.pattern_ID

                blah = np.array([len(xys) for _ in xys])

                self.pattern_dict['pattern_%d' % self.pattern_ID] = {'Rot_mat': rmat.tolist(),
                                                                     'Goodness': float(np.mean(goodness)),
                                                                     'Dist': float(np.mean(dist_loss)),
                                                                     'Spec_Orientation': [s.tolist() for s in spec_ori],
                                                                     'Center_Frame': self.group_dict[group][
                                                                         'Center_Frame'],
                                                                     'Count': conf,
                                                                     'RMS': float(np.mean(rms)),
                                                                     'Num_Peaks': float(np.mean(num_peaks)),
                                                                     'Pos': self.group_dict[group]['Pos'],
                                                                     'xys': xys[np.argmax(blah)]}

                self.pattern_ID += 1
        if not substrate:
            with open(params['peak_dict'], 'w') as json_file:
                json.dump(self.peak_dict, json_file)

            with open(params['pattern_dict'], 'w') as json_file:
                json.dump(self.pattern_dict, json_file)
        else:
            with open(params['peak_dict'], 'w') as json_file:
                json.dump(self.peak_dict, json_file)

            with open(params['pattern_dict'], 'w') as json_file:
                json.dump(self.pattern_dict, json_file)
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
        # pritn(peak_list_r)
        peak_list_f = [xy for xy in zip(xtalsimu.xs, xtalsimu.ys)]

        # print(peak_list_f)
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
        find_rlv = re.compile(r"recip_lattice\d\s+{(.+)}")
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
