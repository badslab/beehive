
import logging
import warnings

import click

logging.getLogger('streamlit.runtime.caching.cache_data_api').setLevel(logging.CRITICAL)

import termite.cli.db
import termite.h5ad

from termite import db


logging.basicConfig(level=logging.INFO)


@click.group("cli")
def cli():
    pass


@cli.command("test")
def test():
    print('?')
    print(db.find_gene_candidates('h.colmg.3', 'raw'))
    

cli.add_command(termite.h5ad.h5ad)
cli.add_command(termite.cli.db.db_group)



def run():
    cli()
