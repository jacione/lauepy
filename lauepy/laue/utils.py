import yaml
from pathlib import Path


def read_config(yml_file):
    with open(yml_file, 'r') as f:
        config = yaml.safe_load(f)
    working_dir = Path(config['working_dir'])
    id_dir = Path(f"{config['working_dir']}/{config['working_id']}")
    if not id_dir.exists():
        id_dir.mkdir(parents=True)
    return config
