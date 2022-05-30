
import logging

import typer

from beehive import config

from beehive.app import gene_expression
from beehive.app import dframes

logging.basicConfig()
app = typer.Typer()



app.add_typer(gene_expression.app, name="gex")
app.add_typer(dframes.app, name="df")


def run():
    app()
