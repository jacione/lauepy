import yaml


def read_config(yml_file):
    with open(yml_file, 'r') as f:
        config = yaml.safe_load(f)
    return config
