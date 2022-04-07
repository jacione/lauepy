"""
A general script for running the AutoLaue package.

Written by Nick Porter
jacioneportier@gmail.com
"""
import click
import lauepy.laue.image_prep as prep
import lauepy.laue.utils as ut
import lauepy.laue.peaks as pk
import lauepy.laue.auto_laue as al
import lauepy.laue.make_grain_dict as grain


@click.command()
@click.argument('config', default='')
def main(config):
    config = ut.read_config(config)
    ut.purge(config)

    # Image prep
    prep.extract_substrate(config)
    prep.cleanup_images(config)

    # Peak indexing
    peak_dict = pk.find_substrate_peaks(config, {})
    peak_dict = pk.find_sample_peaks(config, peak_dict)
    peak_dict = pk.record_positions(config, peak_dict)
    pk.save_peaks(config, peak_dict)

    # Laue analysis
    sim = al.AutoLaue(config)
    sim.index()

    # Grain mapping and spec output
    grain.make_grain_dict(config)


if __name__ == '__main__':
    main()
