import click
import lauepy.laue.auto_laue as al
import lauepy.laue.utils as ut
import lauepy.laue.make_grain_dict as grain


@click.command()
@click.argument('config')
def main(config):
    config = ut.read_config(config)
    sim = al.AutoLaue(config)
    sim.index()
    grain.make_grain_dict(config)


if __name__ == '__main__':
    main()
