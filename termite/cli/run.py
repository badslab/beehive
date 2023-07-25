
import logging

import click

import termite.cli.db
import termite.h5ad
import termite.diffexp
from termite import db

logging.basicConfig(level=logging.INFO)


@click.group("cli")
def cli():
    pass
    

cli.add_command(termite.h5ad.prepare)
cli.add_command(termite.h5ad.h5ad_import)
cli.add_command(termite.diffexp.de_run)
cli.add_command(termite.cli.db.db_group)



def run():
    cli()
