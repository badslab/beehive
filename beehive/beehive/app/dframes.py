"""helpers for the gene expression app."""

import time


import typer
from typer import echo
from typing import List, Optional


from pathlib import Path

from beehive.h5ad_lite import H5adLite

app = typer.Typer()


@app.command()
def h5ad(h5ad_file: Path = typer.Argument(..., exists=True),
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

    study_md = adata.uns['study_md']

    dfx = adata.to_df()
    obs = adata.obs

    study_md['meta'] = {}
    study_md['dimred'] = []

    for k, v in obs.iteritems():

        dtype = 'numerical'

        #polars/ parquest does not like categories
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
        print(obs.shape, oo.shape)

    obs = pd.concat([obs] + obsms, axis=1)

    obs.index.name = '_cell'
    obs = obs.reset_index()

    pl.DataFrame(obs).write_parquet(outbase + '.obs.prq')
    pl.DataFrame(dfx).write_parquet(outbase + '.X.prq')

    with open(outbase + '.yaml','w') as F:
        yaml.dump(study_md, F, Dumper=yaml.SafeDumper)
