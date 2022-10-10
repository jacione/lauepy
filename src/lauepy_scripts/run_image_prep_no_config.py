import click
import src.lauepy.image_prep as prep
from tkinter import Tk
from tkinter.filedialog import askdirectory


@click.command()
@click.option('--data_dir', '-in', help="Path to a directory containing raw images to process")
@click.option('--output-dir', '-out', help="Processed images will be saved in OUTPUT_DIR/clean_images/img_xxxxx")
@click.option('--gaussian-sigma', '-g', default=0.75, help="Sigma to use for Gaussian sigma prior to rolling ball")
@click.option('--rollingball-radius', '-r', type=float, multiple=True,
              help="Rolling ball radius to use for background removal (if given multiple times, RB filters will be "
                   "applied sequentially.)")
def main(data_dir, output_dir, gaussian_sigma, rollingball_radius):
    """

    """
    Tk().withdraw()
    if data_dir is None:
        data_dir = askdirectory()
    if output_dir is None:
        output_dir = askdirectory()
    prep.cleanup_images(data_dir=data_dir, working_dir=output_dir, prep_sample_sigma=gaussian_sigma,
                        prep_sample_radii=rollingball_radius)


if __name__ == '__main__':
    main()
