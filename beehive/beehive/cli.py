

import typer

from beehive import config


app = typer.Typer()

from beehive.app import gene_expression


app.add_typer(gene_expression.app, name="gex")


@app.command()
def help():
    print('help')


def run():
    app()
