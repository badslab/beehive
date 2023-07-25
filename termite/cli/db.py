
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
    print(rv)


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
