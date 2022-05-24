

from functools import partial, lru_cache

import pandas as pd
import polars as pl
import yaml

from beehive import util


DATASETS = {}

diskcache = partial(util.diskcache, where=util.get_datadir('cache'),
                    refresh=True)


@lru_cache(1)
def get_datasets():
    """Return a dict with all dataset."""
    datadir = util.get_datadir('h5ad')
    for yamlfile in datadir.glob('*.yaml'):
        basename = yamlfile.name
        basename = basename.replace('.yaml', '')
        with open(yamlfile, 'r') as F:
            y = yaml.load(F, Loader=yaml.SafeLoader)
        DATASETS[basename] = y

    return DATASETS


def get_dataset(dsid):
    return get_datasets()[dsid]


def get_gene_meta_agg(dsid:str,
                      gene:str,
                      meta:str,
                      nobins:int = 8):
    """Return gene and observation."""
    genedata = get_gene(dsid, gene)
    metadata = get_meta(dsid, meta, nobins=nobins)

    if genedata is None:
        return None
    if metadata is None:
        return None

    rv = (pl.concat([genedata, metadata],
                    how='horizontal')
          .groupby(meta)
          .agg([
              pl.count(),
              pl.mean(gene).alias('mean'),
              pl.std(gene).alias('std'),
              pl.median(gene).alias('median'),
              pl.quantile(gene, 0.25).alias('q25'),
              pl.quantile(gene, 0.75).alias('q75'),
              pl.quantile(gene, 0.01).alias('q01'),
              pl.quantile(gene, 0.99).alias('q99'),
              pl.min(gene).alias('min'),
              pl.max(gene).alias('max'),
          ])
          )

    rv = rv.to_pandas()
    rv = rv.rename(columns={meta: 'cat_value'})
    return rv


def get_gene(dsid, gene):
    """Return expression values for this dataset."""
    datadir = util.get_datadir('h5ad')
    try:
        rv = pl.read_parquet(datadir / f"{dsid}.X.prq", [gene])
    except pl.exceptions.SchemaError:
        return None
    return rv


def get_meta(dsid, col, nobins=8):
    """Return one obs column."""
    ds = get_dataset(dsid)
    dscol = ds['meta'][col]
    datadir = util.get_datadir('h5ad')
    rv = pl.read_parquet(datadir / f"{dsid}.obs.prq", [col])

    if dscol['dtype'] == 'numerical':
        rvq = pd.qcut(rv.to_pandas()[col], nobins, duplicates='drop').astype(str)
        rv[col] = list(rvq)

    return rv


@diskcache()
def get_genes(dsid):
    """Return a list fo genes for this datset."""
    datadir = util.get_datadir('h5ad')
    X = pl.scan_parquet(datadir / f"{dsid}.X.prq")
    return X.columns


@diskcache()
def obslist(dsid):
    """Return a list fo obs columns for this datset."""
    datadir = util.get_datadir('h5ad')
    X = pl.scan_parquet(datadir / f"{dsid}.obs.prq")
    return X.columns
