
import logging

import typer

from beehive import config

from beehive.app import gene_expression
from beehive.app import h5ad
from beehive.app import data

logging.basicConfig(format="%(message)s")
logging.getLogger('beehive').setLevel(logging.INFO)

app = typer.Typer()

app.add_typer(gene_expression.app, name="gex")
app.add_typer(h5ad.app, name="h5ad")
app.add_typer(data.app, name="data")


def run():
    app()
