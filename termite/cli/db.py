
from pathlib import Path    
import os
import click

from termite import db


@click.group("db")
def db_group():
    pass


@db_group.command("head")
@click.argument('table')
def db_head(table):
    head = db.raw_sql(
        f"SELECT * FROM {table} LIMIT 5")
    print(head)


@db_group.command("forget")
@click.argument("experiment")
def forget(experiment: str) -> None:
    "Forget all about one experiment"
    db.forget(experiment)

    
@db_group.command("sql")
@click.argument("sql", nargs=-1)
def db_sql(sql: str) -> None:
    "Run sql."
    sql = ' '.join(sql)
    rv = db.raw_sql(sql)


@db_group.command("status")
def db_status() -> None:
    """Show some stats & table counts."""
    dbfile = os.environ['TERMITE_DB']
    dbsize = os.path.getsize(dbfile)
    print(f"{'db file':<20s} : {dbfile}")
    print(f"{'db size':<20s} : {dbsize:>14_d}")
    tablecount = db.all_table_count()
    for t in sorted(tablecount):
        c = tablecount[t]
        print(f"{t:<20s} : {c:>14_d}")

        
@db_group.command("experiments")
def experiments() -> None:
    """List known experiments."""
    exps = db.get_experiments()
    print(exps)

    
@db_group.command("helptables")
def db_helptables() -> None:
    """Create helper tables."""
    db.helptables()
