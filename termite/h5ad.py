
from functools import lru_cache
import logging
import os
import re
from typing import Dict

import click
import pandas as pd

from termite import db


lg = logging.getLogger(__name__)


@click.group('h5ad')
def h5ad():
    pass


@lru_cache(16)
def get_metadata(h5adfile: str) -> Dict[str, pd.DataFrame]:
    
    rv = {}
    
    obscol_mdfile = h5adfile.replace('.h5ad', '') + '.obscol.tsv'
    if os.path.exists(obscol_mdfile):
        rv['obscol'] = pd.read_csv(obscol_mdfile, sep="\t", index_col=0)
        
    return rv

COMMON_COLUMN_MAPPINGS = """
leiden_scVI | scVI cluster
""".strip().split("\n")
COMMON_COLUMN_MAPPINGS = \
    {a.strip():b.strip() for (a,b)
     in [x.split('|', 1) for x in COMMON_COLUMN_MAPPINGS]
     }

    
def get_obs_column_metadata(h5adfile, adata, column):

    all_md = get_metadata(h5adfile)
    md_obscol = all_md.get('obscol')

    # see if the column is in the obsocl metadata file
    if md_obscol is not None and column in md_obscol['name'].values:
        rv = md_obscol.loc[md_obscol['name'] == column].iloc[0].to_dict()
        return rv
        
    # nothing here? Guess the format from the raw data
    oco = adata.obs[column]
    example=",".join(map(str, oco.head(4)))
    
    if str(oco.dtype) in ['category', 'object']:
        # guess categorical
        uniq = len(oco.unique())
        example = ",".join(map(str, oco.value_counts().sort_values()[:4].index))
        if uniq > 15:
            dtype = 'skip'
        else:
            dtype = 'categorical'
    elif str(oco.dtype).startswith('int'):
        dtype = 'int'
    elif str(oco.dtype).startswith('bool'):
        dtype = 'int'
    elif str(oco.dtype).startswith('float'):
        dtype = 'float'
    else:
        lg.warning(f"Unknown data type for {column} - {oco.dtype}")
        dtype = 'skip'

    if column in COMMON_COLUMN_MAPPINGS:
        nice = COMMON_COLUMN_MAPPINGS[column]
    else:
        nice = " ".join(
            [x.capitalize() for x in re.split(r'[\s_\-]', column)])
        
    return dict(name=column, type=dtype, nice=nice, example=example)

    
@h5ad.command('prepare')
@click.option('-a', '--author', default='unknown')
@click.option('-t', '--title', default='unknown')
@click.option('-s', '--slug')
@click.option('-S', '--study')
@click.option('-v', '--version', type=int, default=1)
@click.option('-o', '--organism', default='unknown')
@click.argument('h5adfile')
def prepare_meta(h5adfile: str,
                 author: str,
                 slug: str,
                 study: str,
                 version: int,
                 organism: str,
                 title: str) -> None:

    import scanpy as sc
    
    basename = h5adfile.replace('.h5ad', '')
    
    if slug is None:
        slug = os.path.basename(h5adfile).replace('.h5ad', '')
        
    if re.match('^[hmx]\.\w{3,10}\.\d{1,5}.*$', slug):
        if organism is 'unknown':
            if slug.startswith('h.'):
                organism = 'human'
            elif slug.startswith('m.'):
                organism = 'mouse'
            elif slug.startswith('x.'):
                organism = 'mixed'

        if study is None:
            study = slug.split('.')[1]

    if study is None:
        study = slug
        
    experiment_md = dict(
        slug = slug,
        author = author,
        study = study,
        version = version,
        title = title,
        organism = organism)

    experiment_md = pd.Series(experiment_md)

    lg.info(f"Loading data: {h5adfile}")
    adata = sc.read_h5ad(h5adfile)

    md_obscol = pd.DataFrame(columns = ['name', 'type', 'nice', 'example'])
    md_obscol = md_obscol.reset_index(drop=True)
    obs = adata.obs

    # regular obs columns
    lg.info("Processing obs table")
    for column in obs.keys():
        if column.startswith('_'):
            continue

        md_obscol.loc[len(md_obscol)] \
            = get_obs_column_metadata(h5adfile, adata, column)

    lg.info("Saving MD files")
    md_obscol.to_csv(f"{basename}.obscol.tsv", sep="\t")

    if not os.path.exists(f"{basename}.experiment.tsv"):
        experiment_md.to_csv(
            f"{basename}.experiment.tsv", sep="\t",
            header=None)

    

@h5ad.command('import')
@click.argument('h5adfile')
@click.option('--expname')
@click.option('--datatype', default='raw', type=str)
@click.option('--layer', default=None, type=str)
def load_import(h5adfile, expname, datatype, layer):

    import scanpy as sc
    import pandas as pd
    conn = db.get_conn()
    
    adata = sc.read_h5ad(h5adfile)

    if expname is None:
        expname = os.path.basename(h5adfile).replace('.h5ad', '')
        
    lg.info(f"data name : {expname}")
    lg.info(f"data gene : {adata.shape[0]:_d}")
    lg.info(f"data cell : {adata.shape[1]:_d}")
    lg.info(f"data type : {datatype}")
    
    total_obs_recs = 0

    
    def _store_col(colname: str,
                   col: pd.DataFrame,
                   dtype: str):
        
        col.columns = ['value']
        
        col.index.name = 'obs'
        col = col.reset_index()
        col['name'] = colname
        col['experiment'] = expname
        
        if dtype == 'skip':
            pass
        elif dtype in ['int', 'float', 'dimred']:

            col['value'] = col['value'].astype(float)
            if db.table_exists('obs_num'):
                conn.sql(f"""
                    DELETE FROM obs_num 
                     WHERE name='{colname}' 
                       AND experiment='{expname}'
                """)
            db.create_or_append('obs_num', col)
        else:
            col['value'] = col['value'].astype(str)

            
            if db.table_exists('obs_cat'):
                conn.sql(f"""
                    DELETE FROM obs_cat 
                     WHERE name='{colname}' 
                       AND experiment='{expname}'
                """)
            db.create_or_append('obs_cat', col)


    for obsmname in adata.obsm_keys():
        lg.info(f"storing dimred {obsmname}")
        if obsmname.startswith('_'):
            continue
        obsmdata = pd.DataFrame(adata.obsm[obsmname][:, :2])
        obsmdata.index = adata.obs_names
        _store_col(f"{obsmname}/0", obsmdata.iloc[:,[0]], 'dimred')
        _store_col(f"{obsmname}/1", obsmdata.iloc[:,[1]], 'dimred')
        
        
    for colname in adata.obs:
        if colname.startswith('_'):
            continue
        col = adata.obs[[colname]].copy()
        colinfo = get_obs_column_metadata(h5adfile, adata, colname)
        if colinfo['type'] == 'skip':
            continue
        
        nice = colinfo['nice']
        dtype = colinfo['type']
        lg.info(f'processing obs | {dtype:14s} | {colname}')
        total_obs_recs += len(col)
        _store_col(nice, col, dtype)

        
    lg.info(f"total obs records {total_obs_recs}")

    
    lg.info("Start storing expression matrix")
    x = adata.to_df(layer=layer)
    x.index.name = 'obs'
    x.columns.name = 'gene'
    melted = x.reset_index().melt(id_vars='obs')
    melted['experiment'] = expname
    melted['datatype'] = datatype
    lg.info(f"melt {melted.shape[0]:_d}")
        
    #remove old data
    lg.info("remove old data")
    if db.table_exists('expr'):
        sql = f"""
            DELETE FROM expr 
             WHERE datatype='{datatype}' 
               AND experiment='{expname}'
            """
        conn.sql(sql)
        
    db.create_or_append('expr', melted)
    
    tablecount = db.all_table_count()
    for t in sorted(tablecount):
        c = tablecount[t]
        print(f"{t:<20s} : {c:>14_d}")

