import re
from pathlib import Path


# inputs
# filename - can include full pathname
# a, b, c lattice params
# alpha, beta, gamma  lattice angles
# inp - in-place vector
# outp - out of place vector


def index_to_macro(index_file, save_dir, name):
    with open(index_file) as f:
        lines = f.readlines()

    lattice_params = None
    recip_lattices = []

    for line in lines:
        if line.startswith('$latticeParameters'):
            lattice_params = re.search('{.*}', line).group(0)
            lattice_params = [float(x) for x in re.split(", ", lattice_params[2:-2])]
            a, b, c, d, e, f = lattice_params
            lattice_params = [a/1e9, b/1e9, c/1e9, d, e, f]
            continue
        if line.startswith('$recip_lattice'):
            recip_lattices.append(re.search("{.*}", line).group(0))
    if lattice_params is None:
        raise IOError("Could not find lattice parameters in Index.txt!")
    for i, text in enumerate(recip_lattices):
        text = re.sub('}{', ",", text[2:-2])
        text = re.split(',', text)
        hklin = [float(text[n]) for n in [0, 3, 6]]
        hklout = [float(text[n]) for n in [1, 4, 7]]
        grain_to_macro(str(Path(save_dir)/f"{name}_{i:03}.mac"), lattice_params, (hklin, hklout))


def grain_to_macro(save_as, lattice_params, orientation):

    a, b, c, alpha, beta, gamma = lattice_params
    a, b, c = a * 1e10, b * 1e10, c * 1e10
    inp, outp = orientation

    with open(save_as, "w") as f:
        f.write("U[\"0\"] = {0:.3f}\n".format(a))  # lattice params
        f.write("U[\"1\"] = {0:.3f}\n".format(b))  # lattice params
        f.write("U[\"2\"] = {0:.3f}\n".format(c))  # lattice params
        f.write("U[\"3\"] = {0:.3f}\n".format(alpha))  # lattice angles
        f.write("U[\"4\"] = {0:.3f}\n".format(beta))  # lattice angles
        f.write("U[\"5\"] = {0:.3f}\n".format(gamma))  # lattice angles

        f.write("U[\"12\"] = {0:.3f}\n".format(inp[0]))  # -8.37  #in-plane HKL
        f.write("U[\"13\"] = {0:.3f}\n".format(inp[1]))  # 0.132
        f.write("U[\"14\"] = {0:.3f}\n".format(inp[2]))  # 7.99
        f.write("U[\"15\"] = {0:.3f}\n".format(outp[0]))  # -0.09  #out-plane HKL
        f.write("U[\"16\"] = {0:.3f}\n".format(outp[1]))  # 11.57
        f.write("U[\"17\"] = {0:.3f}\n".format(outp[2]))  # -0.29

        f.write("U[\"18\"] = 20\n")  # leave all of this as is but actually include them
        f.write("U[\"19\"] = 10\n")  # these define the diffractometer
        f.write("U[\"20\"] = 90\n")  # in/out-plane angles
        f.write("U[\"21\"] = 0\n")
        f.write("U[\"22\"] = 0\n")
        f.write("U[\"23\"] = 0\n")
        f.write("U[\"24\"] = 0\n")
        f.write("U[\"25\"] = 0\n")
        f.write("U[\"26\"] = 90\n")
        f.write("U[\"27\"] = -10\n")
        f.write("U[\"28\"] = 0\n")
        f.write("U[\"29\"] = 20\n")
        f.write("calcG\n")


def write_calibration():
    pass
