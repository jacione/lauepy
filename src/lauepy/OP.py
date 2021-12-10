import numpy as np
import pandas as pd
import scipy.optimize as so
from rxlibs.xmd34 import geometry as geo
from rxlibs.xmd34 import lattice as latt
from scipy.spatial.distance import cdist

import pflibs as pf
import rxlibs_FastForward3 as AP3
import rxlibs_forward_problem_python3_AP as AP

# Set dataset for RP optimization
k = 0

threshold = 20000  # threshold for peak search
tolerance = 1  # baseline for optimization

# define Si crystal
lc_Si = 5.431020511e-10  # lattice constant for Silicon
diamond_pos_list = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0],
                [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25]]
Si_atom_list = [latt.AtomInCell('Si', *atom_pos) for atom_pos in diamond_pos_list]
Si_xtal = latt.Xtal(lc_Si, lc_Si, lc_Si, 90.0, 90.0, 90.0, atomlist=Si_atom_list)
diam_array = np.zeros([len(diamond_pos_list), 3])
for i in range(len(diamond_pos_list)):
    diam_array[i, :] = diamond_pos_list[i]
# rotation and translation vectors (detector)
rot_vec = np.array([1.216432, 1.206319, 1.214645])
trans_vec = np.array([17.317, 0.63, -27.608])
# rot_vec = np.array([1.23420582, 1.20364254, 1.22008472])
# trans_vec = np.array([17.04863603, -0.2922774, -27.68551679])
# rot_vec = np.array([1.23808322, 1.19816954, 1.21511243])
# trans_vec = np.array([16.54574143,  -0.32002094, -27.73272816])
print("Starting R and P")
print('R = ' + str(rot_vec.tolist()))
print('P = ' + str(trans_vec.tolist()))
ccd1 = geo.Detector(516, 517, 0.055, 0.055, 'timepix')
geo_ccd = geo.DetectorGeometry('ccd1', ccd1, trans_vec, rot_vec)

# Import B matrix
B = pd.read_csv('Bmatrix.csv', header=None, delimiter=',').values
# Import experimental Laue patterns and create peak lists
data = pd.read_csv('data2BraggPeakMeasurement.csv', sep=',', index_col=0)

# filename = data['filename'][k] # Filename for optimization
filename = './Staff21-1h_S0029/Staff21-1h_S0029_00000.tif'
[peak_dict, FullPeakList, pre_processed] = pf.get_peaklist(filename, threshold)
# Apply rotation matrix from 2-Bragg peak measurement and overlay simulation with data
# Orientation of the lattice when the diffractometer is at phi = 0, chi = 90, theta = 0
# rotates to the phi, chi, theta of the measured Bragg reflection
delg_gonio = pd.read_csv('delg_gonio' + str(data['H'][k]) + str(data['K'][k]) + str(data['L'][k]) + '.csv', header=None,
                         delimiter=',').values
xtal_rmat = delg_gonio @ B  # rotation matrix
xtal_rmat = np.array(
    [xtal_rmat[i] / np.linalg.norm(xtal_rmat[i]) for i in (0, 1, 2)])  # normalize vectors in rot matrix
energy_cutoff_keV = [10.0, 20.0]  # define energy range
xtal_simu = AP.LaueSimulation(geo_ccd, Si_xtal, xtal_rmat, keV_cutoffs=energy_cutoff_keV, plot_range_padding_percent=0)
#
# fig = plt.figure(figsize=(10, 10), dpi=50)
# ax = fig.add_subplot(111)
# ax.imshow(pre_processed, cmap='terrain_r', vmin=2500, vmax=4500, interpolation='nearest')
xtal_simu.gen_peaks_to_plot()
# xtal_simu.make_plots(ax)


# plt.scatter(FullPeakList['x0'], FullPeakList['y0'], alpha=0.5, edgecolor='yellow', facecolor='yellow', s=81)
# experimentally measured peaks and identified by blob finder

xy_exp = [peak_dict['peak_' + str(k)]['XY'] for k in range(1, len(peak_dict) + 1)]
# simulated peaks by forward problem
xy_sim = [xtal_simu.peaks_plot[0][k].det_coords['ccd1'].tolist() for k in range(0, len(xtal_simu.peaks_plot[0]))]
hkl = [xtal_simu.peaks_plot[0][k].hkl.tolist() for k in range(0, len(xtal_simu.peaks_plot[0]))]
# print("hkl in simulated peaks",hkl)
# print(len(hkl))
q = [xtal_simu.peaks_plot[0][k].q.tolist() for k in range(0, len(xtal_simu.peaks_plot[0]))]
q_norm = [(q[k] / np.linalg.norm(q[k])).tolist() for k in range(0, len(q))]
keV = [xtal_simu.peaks_plot[0][k].keV for k in range(0, len(xtal_simu.peaks_plot[0]))]
# Create DataFrame with hkls, pixel coordinates, and energies for each peak
peaks = [{'H': int(hkl[k][0]), 'K': int(hkl[k][1]), 'L': int(hkl[k][2]), 'px': xy_sim[k][0], 'py': xy_sim[k][1],
          'Qx': q_norm[k][0], 'Qy': q_norm[k][1], 'Qz': q_norm[k][2], 'keV': keV[k]} for k in range(0, len(hkl))]
peaks = pd.DataFrame(data=peaks)
# print(peaks)
pd.DataFrame.to_csv(peaks, 'Peaks' + '.csv', sep=',', header=True, mode='w')
hkl_labels = hkl

# optimization
geoccd = geo.DetectorGeometry('ccd1', ccd1, trans_vec, rot_vec)
xtalsimu = AP3.FastLaueSimulation_list(hkl_labels, geoccd, Si_xtal, xtal_rmat)
peak_list_r = np.array(xy_exp)


def loss_function_distance(x):
    tvec = np.array([x[0], x[1], x[2]])
    rvec = np.array([x[3], x[4], x[5]])
    geoccd = geo.DetectorGeometry('ccd1', ccd1, tvec, rvec)
    xtalsimu = AP3.FastLaueSimulation_list(hkl_labels, geoccd, Si_xtal, xtal_rmat)
    xtalsimu.gen_det_peaks()
    peak_list_r = np.array(xy_exp)
    peak_list_f = [xy for xy in zip(xtalsimu.xs, xtalsimu.ys)]
    # print('peak_list_r is',peak_list_r)
    # print('peak_list_f is',peak_list_f)
    d = cdist(peak_list_r, peak_list_f, metric='euclidean')
    indx_min = np.argmin(d, axis=1)
    distance_sum = 0
    for k4 in range(len(d)):
        # print(d[k4][indx_min[k4]])
        distance_sum += d[k4][indx_min[k4]]

    # print("current R and P is", rvec, tvec)
    # print("average distance is", distance_sum / len(d))

    return distance_sum / len(d)


found = False

for P0 in np.linspace(10, 20, 3):
    if found:
        break
    for P1 in np.linspace(-1, 1, 3):
        if found:
            break
        for P2 in np.linspace(-26, -29, 3):
            if found:
                break
            for R0 in np.linspace(1.15, 1.25, 3):
                if found:
                    break
                for R1 in np.linspace(1.15, 1.25, 3):
                    if found:
                        break
                    for R2 in np.linspace(1.15, 1.25, 3):
                        e = so.minimize(loss_function_distance, [P0, P1, P2, R0, R1, R2], method='Nelder-Mead')
                        print(e.fun)
                        if e.fun < tolerance:
                            P = np.array([e.x[0], e.x[1], e.x[2]])
                            R = np.array([e.x[3], e.x[4], e.x[5]])
                            tolerance = e.fun
                            ebest = e
                            found = True
                            break

# print(ebest)

print("R=", R)
print("P=", P)
RP = [{'R': R[k], 'P': P[k]} for k in range(0, len(R))]
RP = pd.DataFrame(RP)
RP = np.transpose(RP)
pd.DataFrame.to_csv(RP, 'RP.csv', sep=',', header=True, mode='w', index=True)

# # calculate in- and out-of-plane vectors from a*, b*, c* vectors
# M = pf.extract_rlv()
# data = pd.read_csv('/Users/pateras/PycharmProjects/pythonProject/data2BraggPeakMeasurement.csv', sep=',', index_col=0)
# [uvw, hkl] = pf.laue_transform(M, data['phi'][k], data['chi'][k], data['theta'][k])
# print('uvw =', uvw)
# print('hkl =', hkl)
