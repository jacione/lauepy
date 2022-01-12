import click
import lauepy.laue.peaks as pk
import lauepy.laue.utils as ut


@click.command()
@click.argument('config')
def main(config):
    config = ut.read_config(config)
    pk.index_substrate(config)
    pass


if __name__ == '__main__':
    main()
