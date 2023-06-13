
import logging


import click

import termite.cli.db
import termite.h5ad

logging.basicConfig(level=logging.INFO)


@click.group("cli")
def cli():
    pass


cli.add_command(termite.h5ad.h5ad)
cli.add_command(termite.cli.db.db_group)

def run():
    cli()
