
import logging

import typer

logging.basicConfig(format="%(message)s | %(name)s")
logging.getLogger('beehive').setLevel(logging.INFO)


def run():
    from beehive import config
    from beehive.app import data, gene_expression, geneset, h5ad, query

    app = typer.Typer(pretty_exceptions_show_locals=False)

    app.add_typer(gene_expression.app, name="gex")
    app.add_typer(h5ad.app, name="h5ad")
    app.add_typer(data.app, name="data")
    app.add_typer(geneset.app, name="geneset")
    app.add_typer(query.app, name="query")

    app()