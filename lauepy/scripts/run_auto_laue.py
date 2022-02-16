import click
import lauepy.laue.auto_laue as al
import lauepy.laue.utils as ut


@click.command()
@click.argument('config')
def main(config):
    config = ut.read_config(config)
    sim = al.AutoLaue(config)
    sim.index()


if __name__ == '__main__':
    main()
