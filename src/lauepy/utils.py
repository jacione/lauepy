import shutil
from pathlib import Path
from tkinter import Tk, filedialog
import re
import json

import numpy as np
import yaml
from xrayutilities.io import spec

# Get the required keys from the default config file
example = Path(__file__).parents[2] / 'config_example/config.yml'
with open(example, 'r') as F:
    REQUIRED = yaml.safe_load(F)

SPEC_KEYS = ['Piezo_X', 'Piezo_Y', 'Piezo_Z', 'Lab_X', 'Lab_Y', 'Lab_Z']
SPEC_DICT = {
    "px": 'Piezo_X',
    "py": 'Piezo_Y',
    "pz": 'Piezo_Z',
    "labx": 'Lab_X',
    "laby": 'Lab_Y',
    "labz":  'Lab_Z'
}


def new_analysis_cli(entries=None):
    print('New analysis...')
    if entries is None:
        year = input('Year: ')
        exp_id = input('Experiment ID (e.g. LauePUP422): ')
        scan = input('Scan number: ')
    else:
        year, exp_id, scan = entries
    exp_dir = Path(f"/home/beams7/CXDUSER/34idc-work/{year}/{exp_id}")
    if not exp_dir.exists():
        print('Could not find experiment! Check that the ID matches your experiment directory name.')
        return
    alt_id = ''
    work_dir = Path(f'{exp_dir}/Analysis/lauepy_output/scan_{scan:0>4}')
    if work_dir.exists():
        count = 0
        while work_dir.exists():
            alt_id = chr(97 + count)
            work_dir = work_dir.parent / f'scan_{scan:0>4}{alt_id}'
            count += 1
    work_dir.mkdir(parents=True)
    txt = example.read_text()
    txt = re.compile('scan:.*').sub(f'scan: {scan}', txt)
    txt = re.compile('alt_id:.*').sub(f'alt_id: {alt_id}', txt)
    txt = re.compile('year:.*').sub(f'year: {year}', txt)
    txt = re.compile('exp_id:.*').sub(f'exp_id: {exp_id}', txt)
    cfg = work_dir / 'config.yml'
    cfg.touch()
    cfg.write_text(txt)
    print(f'Analysis folder created:\n\t{work_dir}')
    return str(cfg)


def new_analysis(year, exp_id, scan, alt_id, config):
    exp_dir = Path(f"/home/beams7/CXDUSER/34idc-work/{year}/{exp_id}")
    work_dir = Path(f'{exp_dir}/Analysis/lauepy_output/scan_{scan:0>4}{alt_id}')

    # Assertions to make sure the analysis directory is valid
    assert exp_dir.exists()
    assert not work_dir.exists()

    work_dir.mkdir(parents=True)
    cfg = work_dir / 'config.yml'
    cfg.touch()

    for key, val in zip(["year", "exp_id", "scan", "alt_id"], [year, exp_id, scan, alt_id]):
        config[key] = val
    save_config(config, cfg)
    return cfg


def file_prompt():
    root = Tk()
    root.withdraw()
    f = filedialog.askopenfilename()
    return f


def read_config(yml_file):
    """
    Reads a configuration file and returns it as a dictionary.

    :param yml_file: path to a YML file
    :type yml_file: str
    :return: config
    :rtype: dict
    """

    # Load the configuration file
    if yml_file == '':
        yml_file = file_prompt()
    with open(yml_file, 'r') as f:
        cfg = yaml.safe_load(f)

    # Check the config to make sure it has all the right fields
    missing_keys = REQUIRED.keys() - cfg.keys()
    if missing_keys:
        raise KeyError(f'The following required fields were not found in the config file: {missing_keys}')

    # Add the directories to the config dictionary
    year = cfg['year']
    exp_id = cfg['exp_id']
    scan = cfg['scan']
    alt_id = cfg['alt_id']
    if alt_id is None or alt_id == 'None':
        alt_id = ''
    spec_seq = cfg['spec_seq']
    cfg['working_dir'] = f"/home/beams7/CXDUSER/34idc-work/{year}/{exp_id}/Analysis/lauepy_output" \
                         f"/scan_{scan:04}{alt_id}"
    cfg['data_dir'] = f"/home/beams7/CXDUSER/34idc-data/{year}/{exp_id}/AD34idcLaue_{exp_id}{spec_seq}" \
                      f"/{exp_id}{spec_seq}_S{scan:04}"
    cfg['spec_file'] = f"/home/beams7/CXDUSER/34idc-data/{year}/{exp_id}/{exp_id}{spec_seq}.spec"
    cfg['lauepy_dir'] = f'{Path(__file__).parents[1]}'

    # Make sure the data and spec files can be found
    if not Path(cfg['data_dir']).exists():
        raise FileNotFoundError('Could not find the specified DATA DIRECTORY. Check config.')
    if not Path(cfg['spec_file']).exists():
        raise FileNotFoundError('Could not find the specified SPEC FILE. Check config')

    # Load up the calibration data
    cal_path = f'{Path(cfg["working_dir"]).parents[1]}/DetectorCalibration/{cfg["calibration"]}'
    if not cal_path.endswith(".json"):
        cal_path += ".json"
    with open(cal_path) as f:
        cal = json.load(f)
    cfg['det_rotation'] = cal['rot_vec']
    cfg['det_translation'] = cal['trans_vec']
    cfg['det_pixels'] = cal['pixels'][0:2]
    cfg['det_pitch'] = cal['pixels'][2:4]
    cfg['det_name'] = cal['det_name']
    for key, val in cal.items():
        cfg[key] = val

    # Create the working directory if needed
    id_dir = Path(cfg['working_dir'])
    for subdir in ['', 'clean_images', 'peaks', 'substrate', 'grains', 'twins']:
        d = id_dir / subdir
        if not d.exists():
            d.mkdir(parents=True)

    return cfg


def save_config(cfg, yml_file):
    f = Path(yml_file)
    if not f.exists():
        f.parent.mkdir(exist_ok=True)
        f.touch()
        f.write_text(example.read_text())
    text = f.read_text()
    for key, val in cfg.items():
        text = re.compile(f'{key}:.*').sub(f'{key}: {val}', text)
    f.write_text(text)
    return


def read_spec_log(config):
    """
    Reads the positions of a single motor for every frame in a scan.

    :param config: configuration dictionary
    :type config: dict
    :return: array of shape (N,) where N is the number of frames in the scan
    :rtype: ndarray
    """
    # The (scan - 1) indexing is needed because the spec file starts at scan 1 (not zero)
    scan = spec.SPECFile(config['spec_file'])[config['scan'] - 1]
    scan.ReadData()
    axes = [v for k, v in SPEC_DICT.items() if k in scan.command]
    arr = np.empty((scan.data.shape[0], len(axes)))

    for i, ax in enumerate(axes):
        try:
            # Try to find the key in the column headers of the spec table.
            arr[:, i] = np.array(scan.data[ax])
        except ValueError:
            # If the key isn't in the header, then it was held constant. In that case, just use the initial value.
            arr[:, i] = scan.init_motor_pos[f'INIT_MOPO_{ax}']
    return axes, arr


def read_spec_init(config, *keys):
    """
    Reads the initial position of a single motor for a single scan.

    :param config: configuration dictionary
    :type config: dict
    :param keys: motor name
    :type keys: str
    :return: initial motor position
    :rtype: ndarray
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
    id_dir = Path(f"{config['working_dir']}")
    if id_dir.exists():
        shutil.rmtree(id_dir)
        for subdir in ['', 'clean_images', 'peaks', 'substrate', 'grains']:
            d = id_dir / subdir
            if not d.exists():
                d.mkdir(parents=True)
    return


if __name__ == '__main__':
    pass
