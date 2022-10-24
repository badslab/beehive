
import logging

import typer

from beehive import config

from beehive.app import gene_expression
from beehive.app import h5ad
from beehive.app import data
from beehive.app import geneset
from beehive.app import query


logging.basicConfig(format="%(message)s")
logging.getLogger('beehive').setLevel(logging.INFO)

app = typer.Typer(pretty_exceptions_show_locals=False)

app.add_typer(gene_expression.app, name="gex")
app.add_typer(h5ad.app, name="h5ad")
app.add_typer(data.app, name="data")
app.add_typer(geneset.app, name="geneset")
app.add_typer(query.app, name="query")


def run():
    app()
