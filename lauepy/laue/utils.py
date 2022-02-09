import shutil
from pathlib import Path

import numpy as np
import yaml
from xrayutilities.io import spec

# Get the required keys from the default config file
example = Path(__file__).parents[2] / 'config_example/config.yml'
with open(example, 'r') as F:
    REQUIRED = yaml.safe_load(F)


def read_config(yml_file):
    """
    Reads a configuration file and returns it as a dictionary.

    :param yml_file: path to a YML file
    :type yml_file: str
    :return: config
    :rtype: dict
    """
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
    cfg['working_dir'] = f"/home/beams7/CXDUSER/34idc-work/{year}/{exp_id}/Analysis/lauepy_output/scan_{scan:04}"
    cfg['data_dir'] = f"/home/beams7/CXDUSER/34idc-data/{year}/{exp_id}/AD34idcLaue_{exp_id}a/{exp_id}a_S{scan:04}"
    cfg['spec_file'] = f"/home/beams7/CXDUSER/34idc-data/{year}/{exp_id}/{exp_id}a.spec"

    if not Path(cfg['data_dir']).exists():
        raise FileNotFoundError('Could not find the specified DATA DIRECTORY. Check config.')
    if not Path(cfg['spec_file']).exists():
        raise FileNotFoundError('Could not find the specified SPEC FILE. Check config')

    id_dir = Path(cfg['working_dir'])
    for subdir in ['', 'clean_images', 'peaks', 'substrate', 'groups', 'grains']:
        d = id_dir / subdir
        if not d.exists():
            d.mkdir(parents=True)

    return cfg


def read_spec_log(config, *keys):
    """
    Reads the positions of a single motor for every frame in a scan.

    :param config: configuration dictionary
    :type config: dict
    :param keys: motor name
    :type keys: str
    :return: array of shape (N,) where N is the number of frames in the scan
    :rtype: ndarray
    """
    # The (scan - 1) indexing is needed because the spec file starts at scan 1 (not zero)
    scan = spec.SPECFile(config['spec_file'])[config['scan'] - 1]
    scan.ReadData()
    arr = np.empty((scan.data.shape[0], len(keys)))
    for i, key in enumerate(keys):
        try:
            # Try to find the key in the column headers of the spec table.
            arr[:, i] = np.array(scan.data[key])
        except ValueError:
            # If the key isn't in the header, then it was held constant. In that case, just use the initial value.
            arr[:, i] = scan.init_motor_pos[f'INIT_MOPO_{key}']
    return arr


def read_spec_init(config, *keys):
    """
    Reads the initial position of a single motor for a single scan.

    :param config: configuration dictionary
    :type config: dict
    :param keys: motor name
    :type keys: str
    :return: initial motor position
    :rtype: float
    """
    scan = spec.SPECFile(config['spec_file'])[config['scan'] - 1]
    if len(keys) == 1:
        ret = scan.init_motor_pos[f'INIT_MOPO_{keys[0]}']
    else:
        ret = [scan.init_motor_pos[f'INIT_MOPO_{key}'] for key in keys]
    return np.array(ret)


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
