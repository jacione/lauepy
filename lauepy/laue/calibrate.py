import json
from itertools import product
import time

from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tif
import yaml

from lauepy.laue import utils as ut
import lauepy.laue.forward_sim as fsim
import lauepy.laue.spec as sp
from lauepy.laue import image_prep
from lauepy.laue import pflibs
from lauepy.laue.write_specorient import grain_to_spec
from lauepy.rxlibs.xmd34 import geometry as geo
from lauepy.rxlibs.xmd34 import lattice as latt

filename = '/home/beams/CXDUSER/34idc-data/2021/LauePUP921/AD34idcLaue_LauePUP921a/LauePUP921a_S1185/LauePUP921a_S1185_00001.tif'
cryst_params_file = 'lauepy/crystals/Si.json'
det_params_template = 'config_example/detcal.yml'
det_params_output = 'detcal/eiger_params_3_22a.json'
threshold_max = 1e6
threshold_min = 1000
Bmat_file = 'Bmatrix-9-27-21.csv'
delg_gonio_file = 'delg_gonio-09-27-21.csv'
energy_cutoff_keV = [6.0, 18.0]
acceptable_pix_distance = 35
spec_file = '/home/beams/CXDUSER/34idc-data/2021/LauePUP921/LauePUP921a.spec'
scan_number = 1185
spec_macro_out_path = '../../../Experiment/Si100substrate-9-27-21.mac'
trans_vec_step = 0.1
rot_vec_step = .05
# trans_vec_guess = [0.8965888621712237,-19.57934741960132,-60.32659006837108]
# rot_vec_guess = [2.2501011727933413,-0.040611445057554063,2.2080207092648627]
trans_vec_guess = [0.8965888621712237, -19.27934741960132, -60.62659006837108]
rot_vec_guess = [2.2501011727933413, -0.010611445057554063, 2.2080207092648627]
num_steps = 7
# recent opimization at z=4 [ 1.45752795e+00 -2.06882674e+01 -5.98733032e+01  2.24945175e+00 1.76302081e-02  2.22776889e+00]


# Crystal Definition

with open(cryst_params_file) as file:
    cryst_params = json.load(file)

pos_list = cryst_params['pos_list']
latts = cryst_params['lattice_params']
material = cryst_params['material']
print(pos_list)
atom_list = [latt.AtomInCell(material, *atom_pos) for atom_pos in pos_list]
print("list", atom_list)
xtal = latt.Xtal(latts[0], latts[1], latts[2], latts[3], latts[4], latts[5], atomlist=atom_list)

with open(det_params_template) as file:
    det_params = yaml.safe_load(file)

pix_x, pix_y = det_params['det_pixels']
pitch_x, pitch_y = det_params['det_pitch']
det_type = det_params['det_name']

# pix_x,pix_y,pitch_x,pitch_y,det_type = det_params['pixels']

detector = geo.Detector(pix_x, pix_y, pitch_x, pitch_y, det_type)

rot_vec = rot_vec_guess  # det_params['rot_vec']
trans_vec = trans_vec_guess  # det_params['trans_vec']

geo_det = geo.DetectorGeometry(det_params['det_name'], detector, trans_vec, rot_vec)

# Import experimental Laue pattern and create peak list

importedImage = tif.imread(filename).astype('float32')

pre_processed = importedImage[:, :, np.newaxis]
pre_processed[pre_processed > threshold_max] = 0

threshold = threshold_min

[peak_dict, FullPeakList, pre_processed] = image_prep.get_peaklist(pre_processed, threshold)
# Apply rotation matrix from 2-Bragg peak measurement and overlay simulation with data
# use B matrix
B = pd.read_csv(Bmat_file, header=None, delimiter=',').values
delg_gonio = pd.read_csv(delg_gonio_file, header=None, delimiter=',').values
xtal_rmat = delg_gonio @ B

# normalize vectors in rot matrix
for i in (0, 1, 2):
    vec = xtal_rmat[i]
    vec = vec / np.linalg.norm(vec)
    xtal_rmat[i] = vec

xtal_simu = fsim.LaueSimulation(geo_det, xtal, xtal_rmat, keV_cutoffs=energy_cutoff_keV, plot_range_padding_percent=0)

fig = plt.figure(figsize=(18, 10), dpi=50)
ax = fig.add_subplot(111)
ax.imshow(pre_processed, cmap='terrain_r', vmin=0, vmax=5000, interpolation='nearest')

xtal_simu.gen_peaks_to_plot()
xtal_simu.make_plots(ax)

plt.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.5, edgecolor='yellow', facecolor='yellow', s=81)
plt.show()

# experimentally measured peaks and identified by blob finder
xy_exp = [peak_dict['peak_' + str(k)]['XY'] for k in range(2, len(peak_dict) + 1)]
# simulated peaks by forward problem
xy_sim = [xtal_simu.peaks_plot[0][k].det_coords[det_params['det_name']].tolist() for k in
          range(0, len(xtal_simu.peaks_plot[0]))]
hkl = [xtal_simu.peaks_plot[0][k].hkl.tolist() for k in range(0, len(xtal_simu.peaks_plot[0]))]
q = [xtal_simu.peaks_plot[0][k].q.tolist() for k in range(0, len(xtal_simu.peaks_plot[0]))]
q_norm = [(q[k] / np.linalg.norm(q[k])).tolist() for k in range(0, len(q))]
keV = [xtal_simu.peaks_plot[0][k].keV for k in range(0, len(xtal_simu.peaks_plot[0]))]
peaks = [{'H': int(hkl[k][0]), 'K': int(hkl[k][1]), 'L': int(hkl[k][2]), 'px': xy_sim[k][0], 'py': xy_sim[k][1],
          'Qx': q_norm[k][0], 'Qy': q_norm[k][1], 'Qz': q_norm[k][2], 'keV': keV[k]} for k in range(0, len(hkl))]
indices = list()
dict_comparison = dict()

# print(xy_exp)

combs = list(product(peaks, xy_exp))


def loss(a, b):
    return np.linalg.norm(np.array(a) - np.array(b))


xp_list = [((c[0]['H'], c[0]['K'], c[0]['L']), c[1]) for c in combs if
           loss([c[0]['px'], c[0]['py']], c[1]) < acceptable_pix_distance]
print('list of peaks for optimization:', xp_list)

# import rxlibs_FastForward3 as AP3

# # optimization
# tolerance = 1
ind = np.random.randint(len(xp_list), size=10)

# filtered_xy_exp = [xp_list[i][1] for i in ind]
# hkl_labels = [xp_list[i][0] for i in ind]
filtered_xy_exp = [x[1] for x in xp_list]
hkl_labels = [x[0] for x in xp_list]


def loss_function_distance(x):
    tvec = np.array([x[0], x[1], x[2]])
    rvec = np.array([x[3], x[4], x[5]])
    geoccd = geo.DetectorGeometry('ccd1', detector, tvec, rvec)
    xtalsimu = fsim.FastLaueSimulation_list(hkl_labels, geoccd, xtal, xtal_rmat)
    xtalsimu.gen_det_peaks()
    peak_list_r = np.array(filtered_xy_exp)
    peak_list_f = [xy for xy in zip(xtalsimu.xs, xtalsimu.ys)]

    d = cdist(peak_list_r, peak_list_f, metric='euclidean')
    indx_min = np.argmin(d, axis=1)
    distance_sum = 0
    for k4 in range(len(d)):
        distance_sum += d[k4][indx_min[k4]]

    # print("current R and P is", rvec, tvec)
    # print("average distance is", distance_sum / len(d))

    return distance_sum / len(d)


rot_vec_mesh_range = rot_vec_step
trans_vec_mesh_range = trans_vec_step

factor = [1]  # ,0.125,0.1,0.05]
minimum_err = 1000
# trans_vec = [0.8965888621712237,-19.27934741960132,-60.62659006837108]
# rot_vec = [2.2501011727933413,-0.010611445057554063,2.2080207092648627]

min_vals = np.concatenate([trans_vec, rot_vec])
print('Loss for current values: %.4f' % loss_function_distance(min_vals))
print('Current values: ', np.round(min_vals, 3))
# tolerance = [1,0.5,0.2,0.1,0.05,0.001]

for i, f in enumerate(factor):
    start = time.time()
    count = 0
    print('trans_vec range:',
          [np.round((min_vals[i] - trans_vec_step * f, min_vals[i] + trans_vec_step * f), 3) for i in range(3)])
    print('rot_vec range:',
          [np.round((min_vals[i] - rot_vec_step * f, min_vals[i] + rot_vec_step * f), 3) for i in range(3, 6)])
    for P0 in np.linspace(trans_vec[0] - trans_vec_mesh_range * f, trans_vec[0] + trans_vec_mesh_range * f, num_steps):

        for P1 in np.linspace(trans_vec[1] - trans_vec_mesh_range * f, trans_vec[1] + trans_vec_mesh_range * f,
                              num_steps):

            for P2 in np.linspace(trans_vec[2] - trans_vec_mesh_range * f, trans_vec[2] + trans_vec_mesh_range * f,
                                  num_steps):

                for R0 in np.linspace(rot_vec[0] - rot_vec_mesh_range * f, rot_vec[0] + rot_vec_mesh_range * f,
                                      num_steps):

                    for R1 in np.linspace(rot_vec[1] - rot_vec_mesh_range * f, rot_vec[1] + rot_vec_mesh_range * f,
                                          num_steps):

                        for R2 in np.linspace(rot_vec[2] - rot_vec_mesh_range * f, rot_vec[2] + rot_vec_mesh_range * f,
                                              num_steps):

                            e = loss_function_distance([P0, P1, P2, R0, R1,
                                                        R2])  # so.minimize(loss_function_distance, [P0, P1, P2, R0, R1, R2], method='Nelder-Mead', tol=.5)
                            count += 1
                            if e < minimum_err:
                                min_vals = [P0, P1, P2, R0, R1, R2]
                                minimum_err = e
                            if count % 10000 == 0:
                                print('step %s out of %s' % (count, num_steps ** 6))
                                # print('P0 = ', e.x[0], 'P1 = ', e.x[1], 'P2 = ', e.x[2])
                                # print('R0 = ', e.x[3], 'R1 = ', e.x[4], 'R2 = ', e.x[5])
                                print('Minimum Loss: %s Corresponding Vals: %s' % (
                                    np.round(minimum_err, 3), np.float32(np.round(min_vals, 3))))
    end = time.time()
    print('\nTime for optimization loop %d: %.2f seconds\n\n' % (i, end - start))

    rot_vec = min_vals[3:]
    trans_vec = min_vals[:3]

# trans_vec = [0.8965888621712237,-19.27934741960132,-60.62659006837108]
# rot_vec = [2.2501011727933413,-0.010611445057554063,2.2080207092648627]

intensity = [1 for i in range(len(FullPeakList['x0']))]
px = FullPeakList['x0']  # .values
py = FullPeakList['y0']  # .values

xys = []
for i in range(len(px)):
    xys.append([px[i], py[i]])
qs = pflibs.calc_gs(xys, rot_vec, trans_vec, det_params['pixels'])
# print("qs are ",qs)
pflibs.write_mat_params(intensity, qs, cryst_params['material'], cryst_params['lattice_params'],
                        cryst_params['space_group'])

# !./eulerlinux -k 24 -c 179  -f Peaks.txt

goodness, error = pflibs.extract_err()

M = pflibs.extract_rlv()
[delta, _, theta, phi, chi, smot, smot_del, detd, det_name, energy] = sp.parse_spec(spec_file, scan_number)
[uvw, hkl] = pflibs.laue_transform(M, phi, chi, theta)

a, b, c, alpha, beta, gamma = cryst_params['lattice_params']
pflibs.write_spec_orient(spec_macro_out_path, a * 1e10, b * 1e10, c * 1e10, alpha, beta, gamma, uvw, hkl)

geo_det = geo.DetectorGeometry(det_params['det_name'], detector, trans_vec, rot_vec)
importedImage = tif.imread(filename).astype('float32')
pre_processed = importedImage[:, :, np.newaxis]
pre_processed[pre_processed > threshold_max] = 0

threshold = threshold_min

[peak_dict, FullPeakList, pre_processed] = pflibs.get_peaklist(pre_processed, threshold)
# Apply rotation matrix from 2-Bragg peak measurement and overlay simulation with data
# use B matrix
B = pd.read_csv(Bmat_file, header=None, delimiter=',').values
delg_gonio = pd.read_csv(delg_gonio_file, header=None, delimiter=',').values
xtal_rmat = delg_gonio @ B

# normalize vectors in rot matrix
for i in (0, 1, 2):
    vec = xtal_rmat[i]
    vec = vec / np.linalg.norm(vec)
    xtal_rmat[i] = vec

xtal_simu = AP.LaueSimulation(geo_det, xtal, xtal_rmat, keV_cutoffs=[8, 24], plot_range_padding_percent=1)

fig = plt.figure(figsize=(18, 10), dpi=50)
ax = fig.add_subplot(111)
ax.imshow(pre_processed, cmap='terrain_r', vmin=0, vmax=5000, interpolation='nearest')

xtal_simu.gen_peaks_to_plot()
xtal_simu.make_plots(ax)

plt.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.5, edgecolor='yellow', facecolor='yellow', s=81)
plt.show()
