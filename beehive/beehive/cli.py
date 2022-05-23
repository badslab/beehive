

import typer

from beehive import config

app = typer.Typer()

from beehive.app import gene_expression
from beehive.app import dframes


app.add_typer(gene_expression.app, name="gex")
app.add_typer(dframes.app, name="df")


def run():
    app()
