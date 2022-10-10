"""
A general script for running the AutoLaue package.

Written by Nick Porter
jacioneportier@gmail.com
"""
import click
from src.lauepy.write_macro import index_to_macro


@click.command()
@click.argument('index', type=str)
@click.argument('output_dir', type=str)
@click.option('--name', '-n', default='grain', show_default=True, help='Filename for the saved macro(s)')
def main(index, output_dir, name):
    """
    Parses INDEX, which must be a file titled 'Index.txt' output by Euler, and creates a spec macro for each pattern
    found. Saves these macros to OUTPUT_DIR under enumerated titles, starting with '<name>_000.mac'
    """
    if not index.endswith('Index.txt'):
        print("ERROR: Index file must be titled 'Index.txt'")
    index_to_macro(index, output_dir, name)


if __name__ == '__main__':
    main()
