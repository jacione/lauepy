import click
import lauepy.laue.image_prep as prep
import lauepy.laue.utils as ut


@click.command()
@click.argument('config', default='')
def main(config):
    config = ut.read_config(config)
    prep.extract_substrate(config)
    prep.cleanup_images(config)


if __name__ == '__main__':
    main()
