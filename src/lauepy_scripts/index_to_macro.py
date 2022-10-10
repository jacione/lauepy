"""
A general script for running the AutoLaue package.

Written by Nick Porter
jacioneportier@gmail.com
"""
import click
from src.lauepy.write_macro import index_to_macro


@click.command()
@click.argument('index')
@click.argument('output')
@click.option('--name', '-n', default='grain')
def main(index, output, name):
    # Grain mapping and spec output
    index_to_macro(index, output, name)


if __name__ == '__main__':
    main()
