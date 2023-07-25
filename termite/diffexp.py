import pickle
from functools import lru_cache
import logging
import os
from pathlib import Path
import re
from typing import Dict

import click
import duckdb
import pandas as pd
import numpy as np

from termite import db, util, config
from termite.h5ad import get_obs_column_metadata, \
    preprocess_catcol

lg = logging.getLogger(__name__)


def run_de(adata, basefolder, expname, colname, layer, method,
           exp_id, nan_filter=['nan']):

    import scanpy as sc
    import pandas as pd

    cache_dir = Path(basefolder) / 'cache'
    if not cache_dir.exists():
        cache_dir.mkdir()
    cache_file = cache_dir / f"{expname}.{colname}.diffexp.pkl"
        
    # ensure there is no log1p check - responsibility of
    # the user..
    if 'log1p' in adata.uns:
        del adata.uns['log1p']
    
    # remove nan's
    unique = adata.obs[colname].unique()
    
    #ensure lowercase
    nan_filter = [str(x).lower() for x in nan_filter]
    tokeep = []
    for u in unique:
        if str(u).lower() not in nan_filter:
            tokeep.append(u)

    # make sure we have something left to compare
    assert len(tokeep) > 1
    
    if cache_file.exists():
        lg.info("get DE results from cache")
        all_de = pd.read_pickle(cache_file)
        
    else:
        lg.info("running DE")        

        sc.tl.rank_genes_groups(
            adata = adata,
            groupby=colname,
            layer=layer,
            groups = tokeep,
            method=method,
            pts=True)
        rgg = adata.uns['rank_genes_groups']

        # calculate mean expression ..
        raw = pd.DataFrame(data=adata.raw.X,
                           index=adata.raw.obs_names,
                           columns=adata.raw.var_names)
        mex = raw[adata.obs[colname].isin(tokeep)].mean()

        all_de = []
        for i, cv in enumerate(sorted(tokeep)):
            
            print(cv)
            de = pd.DataFrame(dict(
                gene = rgg['names'][cv],
                padj = rgg['pvals_adj'][cv],
                pval = rgg['pvals'][cv],
                lfc = rgg['logfoldchanges'][cv],
            ))
            de['slp'] = -(np.log10(de['pval']).clip(-500, 0))
            de['mex'] = list(mex.reindex(de['gene']))
            del de['pval']
            de['exp_id'] = exp_id
            de['colname'] = colname
            de['colval'] = cv
            all_de.append(de)

        all_de = pd.concat(all_de, axis=0).reset_index(drop=True)
        all_de.to_pickle(cache_file)
        
    return all_de

    

@click.command('de_run')
@click.option('-l', '--layer')
@click.option('-m', '--method', default="wilcoxon")
@click.argument('h5adfile')
@click.argument('columns', nargs=-1)
def de_run(h5adfile, columns, layer, method):

    import scanpy as sc
    import pandas as pd

    basename = h5adfile.replace('.h5ad', '')
    basefolder = Path(h5adfile).parent
    
    # start by reading the experimental metadata
    expdata = pd.read_csv(basename + '.experiment.tsv', sep="\t",
                          index_col=0, header=None).T.loc[1]

    expname = expdata.loc['experiment']
    exp_id = db.get_experiment_id(expname)
    lg.info(f"processing experiment {expname} dbid: {exp_id}")
    
    adata = sc.read_h5ad(h5adfile)
    
    if not columns:
        for colname in adata.obs:
            colinfo = get_obs_column_metadata(h5adfile, adata, colname)
            if colinfo['type'] != 'categorical':
                continue
            print(colname)

    lg.info(f"running DE on experiment {expname}")
    for colname in columns:
        lg.info(f'Processing: {colname}')
        
        # to ensure it is the same as the db by
        # prepping
        if not colname in adata.obs:
            print('column not found', colname, "choose from:")
            for x in adata.obs.columns:
                print(f"'{x}'")
            return
        adata.obs[colname] = preprocess_catcol(adata.obs[colname])
        de = run_de(adata=adata, basefolder=basefolder, exp_id=exp_id,
                    expname=expname, colname=colname, layer=layer,
                    method=method, nan_filter=['nan'])


        print(de.head())
        # ensure older versions of this data are removed
        try:
            db.raw_sql(f'''DELETE FROM diffexp
                            WHERE exp_id = {exp_id}
                              AND colname = '{colname}' ''')
        except duckdb.CatalogException:
            # table does not exist? ignore..
            pass
        
        db.create_or_append('diffexp', de)
        
        
