import click
import src.lauepy.peaks as pk
import src.lauepy.utils as ut


@click.command()
@click.argument('config', default='')
def main(config):
    config = ut.read_config(config)
    peak_dict = pk.find_substrate_peaks(config, {})
    peak_dict = pk.find_sample_peaks(config, peak_dict)
    peak_dict = pk.record_positions(config, peak_dict)
    pk.save_peaks(config, peak_dict)


if __name__ == '__main__':
    main()
