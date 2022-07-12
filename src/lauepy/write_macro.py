import json


# inputs
# filename - can include full pathname
# a, b, c lattice params
# alpha, beta, gamma  lattice angles
# inp - in-place vector
# outp - out of place vector

def grain_to_macro(config):
    with open(f'{config["lauepy_dir"]}/crystals/{config["sample"]}.json') as f:
        sample_params = json.load(f)
    with open(f'{config["lauepy_dir"]}/crystals/{config["substrate"]}.json') as f:
        substrate_params = json.load(f)

    with open(f"{config['working_dir']}/grains/grains.json") as f:
        grain_dict = json.load(f)
    for grain in grain_dict:
        if str(grain) == 'substrate':
            a, b, c, alpha, beta, gamma = substrate_params['lattice_params']
        else:
            a, b, c, alpha, beta, gamma = sample_params['lattice_params']
        a, b, c = a * 1e10, b * 1e10, c * 1e10
        inp, outp = grain_dict[grain]['Spec_Orientation']
        filename = f"{config['working_dir']}/grains/{grain}.mac"

        f = open(filename, "w")
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
        f.close()


def write_calibration():
    pass
