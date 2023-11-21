
import os
from pathlib import Path

import click
import pandas as pd

from termite import db


@click.group("db")
def db_group():
    pass


def f_head(table: str, limit: int = 5):
    return db.raw_sql(
        f"SELECT * FROM {table} LIMIT {limit}")


@db_group.command("head")
@click.argument('table')
@click.option('-t', '--transpose', type=bool, is_flag=True, default=False)
def db_head(table, transpose):
    rv = f_head(table)
    if transpose:
        print(rv.T)
    else:
        print(rv)


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
    print(rv)


@db_group.command("describe")
@click.argument("table")
def db_describe(table: str) -> None:
    "Run sql."
    sql = f"describe {table}"
    rv = db.raw_sql(sql)
    print(rv)



def f_status() -> pd.Series:
    rv = {}
    dbfile = rv['dbfile'] = os.environ['TERMITE_DB']
    dbsize = os.path.getsize(dbfile)
    rv['dbsize'] = f"{dbsize:>14_d}"
    tablecount = db.all_table_count()
    for t in sorted(tablecount):
        c = tablecount[t]
        rv[t] = f"{c:>14_d}"
    return pd.Series(rv)




@db_group.command("status")
def db_status() -> None:
    """Show some stats & table counts."""
    print(f_status())

@db_group.command("experiments")
def experiments() -> None:
    """List known experiments."""
    exps = db.get_experiments()
    for x in exps['experiment']:
        print(x)


@db_group.command("obscol")
@click.argument('experiment')
def obscol(experiment: str) -> None:
    """List known experiments."""
    exp_id = db.get_experiment_id(experiment)
    catcol = db.raw_sql(f'''
         SELECT name FROM help_obs_cat
          WHERE exp_id = '{exp_id}'
         ''')['name']
    numcol = db.raw_sql(f'''
         SELECT name FROM help_obs_num
          WHERE exp_id = '{exp_id}'
         ''')['name']
    print(catcol)
    for c in catcol:
        print(f"cat\t{c}")
    for n in numcol:
        print(f"num\t{n}")


@db_group.command("helptables")
def db_helptables() -> None:
    """Create helper tables."""
    db.helptables()
