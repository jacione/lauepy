import shutil
from pathlib import Path
import re
from xrayutilities.io import spec

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


def read_spec_log(config, key):
    scan = spec.SPECFile(config['spec_file'])[config['scan'] - 1]
    scan.ReadData()
    try:
        arr = scan.data[key]
    except KeyError:
        arr = np.ones(scan.data.shape[0]) * scan.init_motor_pos[f'INIT_MOPO_{key}']
    return arr


def read_spec_init(config, key):
    scan = spec.SPECFile(config['spec_file'])[config['scan'] - 1]
    return scan.init_motor_pos[f'INIT_MOPO_{key}']


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
