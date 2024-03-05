"""
Command line interface dispatcher.

"""

# late imports speed up startup time!
# pylint: disable=import-outside-toplevel

from functools import partial
import importlib
import logging
from pprint import pprint   # noqa: F401
import sys
from typing import Any, Union

import click
from click.core import Context
from rich.logging import RichHandler

from . import db, util
from . import metadata_tools as mdtools
from . import cli_query
from . import cli_db_upload


FORMAT = "%(message)s"
if '-vv' in sys.argv:
    loglevel = logging.DEBUG
elif '-v':
    loglevel = logging.INFO
else:
    loglevel = logging.WARNING


logging.basicConfig(
    level=loglevel, format=FORMAT,
    datefmt="[%X]", handlers=[RichHandler()])

lg = logging.getLogger(__name__)
lg.debug("Start CellHive!")


@click.group()
@click.option('-v', 'verbose', count=True)
@click.option('--db', 'dbfile', type=click.Path(exists=False))
@click.pass_context
def cli(ctx: Context, dbfile: str, verbose: int) -> None:
    """Click group."""
    ctx.ensure_object(dict)
    ctx.obj['chdb'] = db.CHDB(dbfile)


@cli.group("db")
@click.pass_context
def db_group(ctx: Context) -> None:
    """Database functions."""


@db_group.command()
@click.pass_context
def status(ctx: Context) -> None:
    """Show db status."""
    chdb = ctx.obj['chdb']
    for k, v in chdb.status().items():
        print(f"{k:<14}: {v:>15_d}")


@db_group.command("sql")
@click.argument("sql", nargs=-1)
@click.option("-T", "transpose",
              is_flag=True, show_default=False, default=False,
              help="transpose output.")
@click.pass_context
def db_sql(ctx: Context, transpose: bool, sql: str) -> None:

    "Run sql."
    sql = ' '.join(sql)
    rv = ctx.obj['chdb'].sql(sql)
    if transpose:
        print(rv.T)
    else:
        print(rv)


@cli.command()
def version() -> None:
    """Print version to screen."""
    print(importlib.metadata.version("cellhive"))


@db_group.command()
@click.pass_context
def meta_tables(
        ctx: Context,
) -> None:

    """(re-)build meta data tables."""

    # database object.
    chdb = ctx.obj['chdb']

    lg.info("Create experiment help table")
    chdb.sql("DROP TABLE IF EXISTS dataset_meta")
    chdb.sql("""
        CREATE TABLE dataset_meta AS
          SELECT dataset_id,
                 count(*) as no_datapoints,
                 count(distinct obs) as no_cells,
                 count(distinct gene) as no_genes
            FROM expr
           GROUP BY dataset_id""")

    lg.info("Create help_gene table")
    chdb.sql("DROP TABLE IF EXISTS gene_meta")
    chdb.sql("""
       CREATE TABLE gene_meta AS
           SELECT distinct dataset_id, gene,
                  SUM(value) AS sumval,
                  SUM(LEAST(value, 1)) / COUNT(value) AS fracnonzero
             FROM expr
            GROUP BY dataset_id, gene
            ORDER BY dataset_id ASC, sumval DESC """)


#add external commands
cli.add_command(cli_query.query)
cli.add_command(cli_db_upload.upload)

def main() -> None:
    """Run the main click CLI function."""
    cli()
