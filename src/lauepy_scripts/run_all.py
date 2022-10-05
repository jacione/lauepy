"""
A general script for running the AutoLaue package.

Written by Nick Porter
jacioneportier@gmail.com
"""
import click
import src.lauepy.image_prep as prep
import src.lauepy.utils as ut
import src.lauepy.peaks as pk
import src.lauepy.auto_laue as al
import src.lauepy.make_grain_dict as grain


@click.command()
@click.argument('config', default='')
def main(config):
    config = ut.read_config(config)
    ut.purge(config)

    # Image prep
    prep.extract_substrate(**config)
    prep.cleanup_images(**config)

    # Peak indexing
    peak_dict = pk.find_substrate_peaks({}, **config)
    peak_dict = pk.find_sample_peaks(peak_dict, **config)
    peak_dict = pk.record_positions(peak_dict, **config)
    pk.save_peaks(peak_dict, **config)

    # Laue analysis
    sim = al.AutoLaue(config)
    sim.index()

    # Grain mapping and spec output
    grain.make_grain_dict(config)


if __name__ == '__main__':
    main()
