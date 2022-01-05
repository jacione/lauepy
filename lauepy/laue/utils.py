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
