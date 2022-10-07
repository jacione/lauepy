from pathlib import Path

import click
import src.lauepy.auto_laue as al
import src.lauepy.utils as ut
import src.lauepy.make_grain_dict as grain
import src.lauepy.find_twins as twins


@click.command()
@click.argument('config', default='')
def main(config):
    config = ut.read_config(config)
    p = Path(f"{config['working_dir']}/grains")
    if 'n' not in input('Run pattern indexing? (Y/n)').lower():
        for f in p.iterdir():
            f.unlink()
        sim = al.AutoLaue(config)
        sim.index()
    grain.make_grain_dict(**config)
    if 'n' not in input('Run twin finding? (Y/n)').lower():
        if twins.find_possible_twins(**config):
            twins.find_twins(**config)
            twins.cleanup_directory()


if __name__ == '__main__':
    main()
