
import click
from click.core import Context


@click.group('q')
@click.pass_context
def query(ctx: Context):
    """Query the database."""


@query.command
@click.pass_context
@click.argument('gene')
def gene(ctx: Context, gene: str):
    chdb = ctx.obj['chdb']

    sql = f"""
        SELECT experiment_md.dataset, gene_meta.fracnonzero
          FROM gene_meta
          JOIN experiment_md ON gene_meta.dataset_id = experiment_md.dataset_id
         WHERE gene = '{gene}'
         ORDER BY fracnonzero DESC
         LIMIT 5
    """
    result = chdb.sql(sql)
    print(result.to_string(index=False))
