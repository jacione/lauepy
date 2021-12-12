import click
from lauepy.laue.background_subtraction import remove_bkg
import lauepy.laue.utils as ut


@click.command()
@click.argument('config')
def main(config):
    config = ut.read_config(config)
    remove_bkg(config)


if __name__ == '__main__':
    main()
