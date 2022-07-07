"""helpers for the gene expression app."""

import logging
import os
import time
from pathlib import Path
from pprint import pprint

import typer
from typer import echo
from typing import List, Optional

from beehive import util, expset

app = typer.Typer()

lg = logging.getLogger(__name__)


@app.command("metadata")
def h5ad_uns(h5ad_file: Path = typer.Argument(..., exists=True)):

    import scanpy as sc
    lg.warning(f"Loading on {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file, backed='r')
    pprint(adata.uns['study_md'])


@app.command("set")
def h5ad_set(h5ad_file: Path = typer.Argument(..., exists=True),
             key: str = typer.Argument(...),
             val: str = typer.Argument(...)):

    import scanpy as sc
    lg.warning(f"Setting on {h5ad_file}")
    lg.warning(f'Changing "{key}" to "{val}" in `.uns["study_md"]`')
    adata = sc.read_h5ad(h5ad_file, backed='r+')
    adata.uns['study_md'][key] = val
    adata.write()


@app.command("del")
def h5ad_del(h5ad_file: Path = typer.Argument(..., exists=True),
             key: str = typer.Argument(...)):

    import scanpy as sc
    lg.warning(f"Changing on {h5ad_file}")
    lg.warning(f'Removing "{key}" from `.uns["study_md"]`')
    adata = sc.read_h5ad(h5ad_file, backed='r+')
    del adata.uns['study_md'][key]
    adata.write()


@ app.command("prepare")
def h5ad_convert(h5ad_file: Path = typer.Argument(..., exists=True),
                 author: str = None,
                 title: str = None,
                 ):
    """Convert to polars/parquet dataframes."""

    import scanpy as sc
    import polars as pl
    import pandas as pd
    import yaml

    outbase = h5ad_file
    lg.info(f"Filename for IO: {outbase}")

    adata = sc.read_h5ad(h5ad_file)

    # automatically fails if not exist
    study_md = adata.uns['study_md']

    for field in ['author', 'title', 'study', 'organism', 'datatype',
                  'short_title', 'year']:
        assert field in study_md

    study_md['year'] = int(study_md['year'])

    try:
        adata.raw.to_adata()
        lg.warning("Adata has a `.raw`! Are you sure you have the correct")
        lg.warning("data in .X?")
    except:  # NOQA: E722
        pass

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
        if (k.endswith('__lfc') or k.endswith('__padj')):
            lg.info(f"Processing var {k}")
            keep.append(k)

    var = var[keep].copy()

    for k, v in var.iteritems():
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

        # polars/parquest does not like categories
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

    lg.info("Writing output files to:")
    lg.info(" - " + str(outbase.with_suffix('.obs.prq')))
    pl.DataFrame(obs).write_parquet(outbase.with_suffix('.obs.prq'))
    lg.info(" - " + str(outbase.with_suffix('.var.prq')))
    pl.DataFrame(var).write_parquet(outbase.with_suffix('.var.prq'))
    lg.info(" - " + str(outbase.with_suffix('.X.prq')))
    pl.DataFrame(dfx).write_parquet(outbase.with_suffix('.X.prq'))

    lg.info(" - " + str(outbase.with_suffix('.yaml')))
    with open(outbase.with_suffix('.yaml'), 'w') as F:
        yaml.dump(study_md, F, Dumper=yaml.SafeDumper)
