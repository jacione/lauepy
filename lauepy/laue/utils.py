import shutil
from pathlib import Path
import re

import yaml
import numpy as np

# Get the required keys from the default config file
example = Path(__file__).parents[2] / 'config_example/config.yml'
with open(example, 'r') as F:
    REQUIRED = yaml.safe_load(F)

FIND_ANGLES = re.compile(r"#P0 (.*?)\n", re.DOTALL)
FIND_PIEZO = re.compile(r"(-?\d+.*?)\n")


def read_config(yml_file):
    # Load the configuration file
    with open(yml_file, 'r') as f:
        cfg = yaml.safe_load(f)

    # Check the config to make sure it has all the right fields
    missing_keys = REQUIRED.keys() - cfg.keys()
    if missing_keys:
        raise KeyError(f'The following required fields were not found in the config file: {missing_keys}')

    # Create the working directory if needed
    year = cfg['year']
    exp_id = cfg['exp_id']
    scan = cfg['scan']
    cfg['working_dir'] = f"/home/beams/CXDUSER/34idc-work/{year}/{exp_id}/Analysis/lauepy_output/scan_{scan:04}"
    cfg['data_dir'] = f"/home/beams/CXDUSER/34idc-data/{year}/{exp_id}/AD34idcLaue_{exp_id}a/{exp_id}a_S{scan:04}"
    cfg['spec_file'] = f"/home/beams/CXDUSER/34idc-data/{year}/{exp_id}/{exp_id}a.spec"

    if not Path(cfg['data_dir']).exists():
        raise FileNotFoundError('Could not find the specified DATA DIRECTORY. Check config.')
    if not Path(cfg['spec_file']).exists():
        raise FileNotFoundError('Could not find the specified SPEC FILE. Check config')

    id_dir = Path(cfg['working_dir'])
    if not id_dir.exists():
        id_dir.mkdir(parents=True)

    return cfg


def get_spec(config):
    with open(config['spec_file']) as f:
        spec_file = f.read()
    scan = config['scan']

    find_scan = re.compile(r"#S %s (.*?)#S %s" % (scan, scan+1), re.DOTALL)
    spec_txt = find_scan.findall(spec_file)
    if len(spec_txt) == 0:
        find_scan = re.compile(r"#S %s .*" % scan, re.DOTALL)
        spec_txt = find_scan.findall(spec_file)

    return spec_txt[0]


def get_spec_angles(config):
    spec_txt = get_spec(config)
    vals = FIND_ANGLES.findall(spec_txt)[0].split(" ")
    theta, chi, phi = np.array(vals[1:4], dtype=float)

    return phi, chi, theta


def get_spec_coords(config):
    spec_txt = get_spec(config)
    spec_txt = spec_txt.split('Monitor  Detector\n')[1].split('#')[0]
    coords = np.array([frame.split(" ")[:2] for frame in FIND_PIEZO.findall(spec_txt)], dtype='f')

    return coords


def purge(config):
    """
    Cleans out the working directory
    """
    id_dir = Path(f"{config['working_dir']}/{config['working_id']}")
    if id_dir.exists():
        shutil.rmtree(id_dir)
        id_dir.mkdir()
    return


if __name__ == '__main__':
    pass
