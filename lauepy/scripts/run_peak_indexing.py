import click
import lauepy.laue.peaks as pk
import lauepy.laue.utils as ut


@click.command()
@click.argument('config')
def main(config):
    config = ut.read_config(config)
    peak_dict = pk.init_peaks()
    peak_dict = pk.find_substrate_peaks(config, peak_dict)
    peak_dict = pk.find_sample_peaks(config, peak_dict)
    peak_dict = pk.record_positions(config, peak_dict)
    pk.save_peaks(config, peak_dict)


if __name__ == '__main__':
    main()
