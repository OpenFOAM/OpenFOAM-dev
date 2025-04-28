# mymodule/__init__.py
# Empty file to make the directory a package

# mymodule/cli.py
import click
from . import melt_front

@click.group()
def cli():
    """My module with nested commands."""
    pass

@cli.command('melt-front')
@click.argument('case_dir')
@click.argument('threshold', default=0.5)
def melt_front_command(case_dir, threshold):
    """Calculate melt front properties and save in postProcessing/"""
    melt_front.main(case_dir, threshold)
