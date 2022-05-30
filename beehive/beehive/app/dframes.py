"""helpers for the gene expression app."""

import logging
import time


import typer
from typer import echo
from typing import List, Optional


from pathlib import Path

from beehive import util, expset

app = typer.Typer()
lg = logging.getLogger(__name__)


@app.command()
def test():
    expset.get_datasets()


@app.command()
def h5ad(h5ad_file: Path = typer.Argument(..., exists=True),
         author: str = None,
         title: str = None,
          ):
    """Convert to polars/parquet dataframes."""

    import scanpy as sc
    import polars as pl
    import pandas as pd
    import yaml

    outbase = h5ad_file
    if str(outbase).endswith('.h5ad'):
        outbase = str(h5ad_file).replace('.h5ad', '')

    adata = sc.read_h5ad(h5ad_file)

    if 'study_md' in adata.uns:
        study_md = adata.uns['study_md']
    else:
        study_md = {}

    if author is not None:
        study_md['author'] = author

    if title is not None:
        study_md['title'] = title

    # try:
    #     adata = adata.raw.to_adata()
    # except:
    #     print("could not switch to raw")

    dfx = adata.to_df()

    obs = adata.obs
    var = adata.var
    var.index = adata.var_names
    var.index.name = 'gene'

    study_md['meta'] = {}
    study_md['diffexp'] = {}
    study_md['dimred'] = []

    keep = []
    for k in var:
        print(k)
        if (k.endswith('__lfc') or k.endswith('__padj')):
            keep.append(k)

    var = var[keep].copy()

    for k, v in var.iteritems():
        print("processing DE", k)
        if not k.endswith('__padj'):
            continue

        kgroup, kkey, _ = k.rsplit('__', 2)
        if kgroup not in study_md['diffexp']:
            study_md['diffexp'][kgroup] = dict(keys=[kkey])
        else:
            study_md['diffexp'][kgroup]['keys'].append(kkey)

        vv = v.sort_values().head(8)
        topgenes = list(vv.index)
        study_md['diffexp'][kgroup]['topgenes'] = topgenes


    for k, v in obs.iteritems():
        dtype = 'numerical'

        #polars/parquest does not like categories
        if str(v.dtype) == 'category':
            obs[k] = v.astype('str')
            dtype = 'categorical'

        study_md['meta'][k] = dict(dtype=dtype)

    obsms = []
    for k in adata.obsm_keys():
        if k.startswith('_'):
            continue
        study_md['dimred'].append(k)

        oo = pd.DataFrame(adata.obsm[k], index=obs.index)
        oo.columns = '_' + k + '_' + oo.columns.astype(str)
        obsms.append(oo)

    obs = pd.concat([obs] + obsms, axis=1)

    obs.index.name = '_cell'
    obs = obs.reset_index()

    var = var.T
    var.index.name = 'field'
    var = var.reset_index()

    pl.DataFrame(obs).write_parquet(outbase + '.obs.prq')
    pl.DataFrame(var).write_parquet(outbase + '.var.prq')
    pl.DataFrame(dfx).write_parquet(outbase + '.X.prq')

    with open(outbase + '.yaml','w') as F:
        yaml.dump(study_md, F, Dumper=yaml.SafeDumper)
