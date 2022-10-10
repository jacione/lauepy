import click
import src.lauepy.image_prep as prep


@click.command()
@click.argument('data_dir')
@click.argument('output_dir')
@click.option('--gaussian-sigma', '-g', default=0.75,
              help="Sigma to use for Gaussian filter prior to rolling ball.")
@click.option('--rollingball-radius', '-r', type=float, multiple=True, default=(45, 15, 5),
              help="Rolling ball radius to use for background removal. (Can be provided multiple times, in which case "
                   "RB filters will be applied sequentially.)")
def main(data_dir, output_dir, gaussian_sigma, rollingball_radius):
    """
    This script applies the LauePy image-processing recipe to all images in DATA_DIR, then saves them in <OUTPUT_DIR>/
    clean_images/img_XXXXX.tiff
    """
    prep.cleanup_images(data_dir=data_dir, working_dir=output_dir, prep_sample_sigma=gaussian_sigma,
                        prep_sample_radii=rollingball_radius)


if __name__ == '__main__':
    main()
