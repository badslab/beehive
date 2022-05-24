"""helpers for the gene expression app."""

import typer
from typer import echo
from typing import List, Optional

from pathlib import Path


app = typer.Typer()


@app.command()
def obs(h5ad_file: Path = typer.Argument(..., exists=True),
        skip: Optional[List[str]] = typer.Option(None),
        ):
    """Show obs fields."""
    import scanpy as sc
    import re

    skipres = []
    for s in skip:
        skipres.append(re.compile(s, re.IGNORECASE))

    adata = sc.read_h5ad(h5ad_file, backed='r')

    for name, dtype in adata.obs.dtypes.iteritems():
        skipthis = False
        for sre in skipres:
            if sre.search(name):
                skipthis = True
                break
        if skipthis:
            continue

        print(name, dtype, dtype.kind)


@app.command()
def load(h5ad_file: Path = typer.Argument(..., exists=True),
         sqlite_db: Path = typer.Argument(...),
         dsid: str = None,
         author: str = None,
         title: str = None,
         skip: Optional[List[str]] = None,
         field: Optional[List[str]] = typer.Option(None),
         categ: Optional[List[str]] = typer.Option(None),
         ):
    """Load h5ad into a sqlite database."""
    import re
    import os
    import sys
    import sqlite3

    import numpy as np
    import pandas as pd
    import scanpy as sc

    if dsid:
        dataset_id = dsid
    else:
        dataset_id = os.path.basename(h5ad_file).replace('.h5ad', '')

    echo(f"loading   : {h5ad_file}")
    echo(f"basename  : {dataset_id}")
    echo("Force cat : " + " ".join(categ))

    skipres = []
    for s in skip:
        echo("skip re: {s}")
        skipres.append(re.compile(s, re.IGNORECASE))

    adata = sc.read_h5ad(h5ad_file)
    sqldb = sqlite3.connect(sqlite_db)

    study_data_raw = {}
    study_data_raw['dataset_id'] = dataset_id

    if 'study_md' in adata.uns:
        aus = adata.uns['study_md']
        study_data_raw['title'] = aus['title']
        study_data_raw['author'] = aus['author']

    if author:
        study_data_raw['author'] = author
    if title:
        study_data_raw['title'] = title

    assert 'author' in study_data_raw
    assert 'title' in study_data_raw

    study_data = pd.DataFrame(pd.Series(study_data_raw)).T

    cursor = sqldb.cursor()
    try:
        cursor.execute(f"DELETE FROM study WHERE dataset_id='{dataset_id}'")
    except sqlite3.OperationalError as e:
        if str(e) != "no such table: study":
            raise

    study_data.to_sql('study', sqldb, if_exists='append', index=False)

    if adata.raw is None:
        af = adata.to_df()
    else:
        af = adata.raw.to_adata().to_df()

    print(f"Data loaded, {af.shape[0]} cells, {af.shape[1]} genes")
    assert list(adata.obs.index) == list(af.index)

    for k in adata.obs:

        if field and k not in field:
            echo(f"Skipping: {k}")
            # if fields are specified - process only this field
            continue

        if k.startswith('_'):
            echo(f"Skipping: {k}")
            continue

        skip = False
        for sre in skipres:
            if sre.search(k):
                skip = True
                break
        if skip:
            echo(f"Skipping - re: {k}")
            continue

        try:
            cursor.execute(
                f"""DELETE FROM data
                     WHERE dataset_id='{dataset_id}'
                       AND cat_name='{k}'
                """)
        except sqlite3.OperationalError as e:
            if str(e) != "no such table: data":
                raise


        bins = None

        if k in categ:
            adata.obs[k] = adata.obs[k].astype(str).astype('category')

        elif str(adata.obs[k].dtype) != 'category':
            if adata.obs[k].dtype.kind not in 'biuf':
                # not categorical, not numerical??
                # ignore
                echo(f"Unknown datatype {k}")
                continue
            # column is numerical - create bins
            nobins = int(len(adata.obs) / 1000)
            nobins = max(3, min(nobins, 10))
            bins = pd.qcut(adata.obs[k], q=nobins,
                           duplicates="drop")

        print("processing", k, end=' : ')

        if bins is not None:
            print('{B}', end=' ')
            grouper = af.groupby(bins)
        else:
            grouper = af.groupby(adata.obs[k])

        collstats = {}
        aggfuncs = dict(mean=np.mean,
                        std=np.std,
                        min=np.min,
                        max=np.max)

        aggfuncs.update({
            f"q{q:02d}": (np.quantile, dict(q=q/100, axis=0))
            for q in [1, 25, 50, 75, 99]})

        for n, func in aggfuncs.items():

            print(n, end=' ')
            sys.stdout.flush()

            if isinstance(func, tuple):
                func, kwargs = func
                agg = grouper.agg(func, **kwargs)
            else:
                agg = grouper.agg(func)

            agg.index.name = 'cat_value'
            agg.columns.name = 'gene'
            agg = agg.reset_index()\
                     .melt(id_vars='cat_value')
            agg = agg.set_index(['cat_value', 'gene'])
            collstats[n] = agg['value']

        print(' saving')
        agg = pd.DataFrame.from_dict(collstats).reset_index()
        agg['cat_name'] = k
        agg['dataset_id'] = dataset_id
        agg['cat_value'] = agg['cat_value'].astype(str)
        agg.to_sql('data', sqldb, if_exists='append', index=False)

    try:
        cursor = sqldb.cursor()

        cursor.execute('''
            CREATE INDEX IF NOT EXISTS "idx_dataset_id"
                 ON data("dataset_id")''')

        cursor.execute('''
            CREATE INDEX IF NOT EXISTS "idx_dataset_cat"
                ON data("dataset_id", "cat_name", "gene")''')

    except sqlite3.OperationalError as e:
        if "no such table" not in str(e):
            # maybe nothing was inserted?
            raise
