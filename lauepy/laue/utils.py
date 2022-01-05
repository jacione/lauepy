import shutil
from pathlib import Path

import yaml

# Get the required keys from the default config file
example = Path(__file__).parents[1] / 'config_example/config.yml'
with open(example, 'r') as f:
    REQUIRED = yaml.safe_load(f)


def read_config(yml_file):
    # Load the configuration file
    with open(yml_file, 'r') as f:
        config = yaml.safe_load(f)

    # Check the config to make sure it has all the right fields
    if missing_keys := REQUIRED.keys() - config.keys():
        raise KeyError(f'The following required fields were not found in the config file: {missing_keys}')

    # Create the working directory if needed
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
