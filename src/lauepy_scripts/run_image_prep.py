import click
import src.lauepy.image_prep as prep
import src.lauepy.utils as ut


@click.command()
@click.argument('config', default='')
def main(config):
    config = ut.read_config(config)
    prep.extract_substrate(**config)
    prep.cleanup_images(**config)


if __name__ == '__main__':
    main()
