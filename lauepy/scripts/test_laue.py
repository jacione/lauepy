import lauepy.laue.forward_sim as fs
from lauepy.rxlibs.xmd34 import lattice as lat
from lauepy.rxlibs.xmd34 import geometry as geo

import numpy as np
from matplotlib import pyplot as plt

xtal_params = {
    "lattice_params": [3.905268e-10, 3.905268e-10, 3.905268e-10, 90, 90, 90],
    "space_group": 221,
    "material": "SrTiO3",
    "pos_list": [["Sr", [0, 0, 0]],
                 ["Ti", [0.5, 0.5, 0.5]],
                 ["O", [0.5, 0.5, 0]]]
}
ccd1 = geo.Detector(1024, 512, 0.075, 0.075, 'eiger')
geo_ccd = geo.DetectorGeometry('ccd1', ccd1,
                               [-0.6040795238420961, -26.651050662120888, -68.15917895706454],
                               [2.238706080280818, -0.0042931382928194685, 2.2207979150987223]
                               )
# rmat = [[0.2673115, -0.6157260, 0.7412328],
#         [-0.5891335, 0.5042874, 0.6313604],
#         [-0.7625393, -0.6054550, -0.2279430]]
# rmat = np.array(rmat).T
# hkls = [[-5, -1, 3], [-7, -1, 3], [-5, 1, 3], [-3, 1, 3], [-6, 0, 4], [-7, 1, 3]]
rmat = [[-0.2101211, -0.3156216, -0.9253281],
        [-0.6588140, -0.6535908, 0.3725360],
        [-0.7223664, 0.6878969, -0.0706027]]
rmat = np.array(rmat).T
hkls = [[1, 0, 1], [3, 0, 2], [3, 0, 4], [5, -1, 5], [3, 0, 5], [10, -2, 7], [10, -1, 7]]
atom_list = [lat.AtomInCell(atom_pos[0], *atom_pos[1]) for atom_pos in xtal_params['pos_list']]
xtal = lat.Xtal(*xtal_params['lattice_params'], atomlist=atom_list)
sim = fs.FastLaueSimulation_list(hkls, geo_ccd, xtal, rmat)
sim.gen_det_peaks()
coords = np.array([sim.xs, sim.ys]).T

sim2 = fs.LaueSimulation(geo_ccd, xtal, rmat, keV_cutoffs=(3, 30))
sim2.gen_peaks_to_plot()
coords2 = np.array([sim2.xs, sim2.ys]).T

for c in coords:
    print(c in coords2)

plt.scatter(sim2.xs, sim2.ys, c='k', label='LaueSimulation')
plt.scatter(sim.xs, sim.ys, edgecolor='red', facecolor='None', s=160, label='FastLaueSimulation')
plt.legend()
plt.show()
