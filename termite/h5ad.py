
import re
import os
import logging

from termite import db

import click

lg = logging.getLogger(__name__)


@click.group('h5ad')
def h5ad():
    pass

@h5ad.group('meta')
def meta():
    pass


def get_obs_column_metadata(adata, column):
    if 'termite_obs_meta' in adata.uns:
        # see if the column is in the adata file
        bom = adata.uns['termite_obs_meta']

        if column in bom['name']:
            return bom.loc[bom['name'] == column].to_dict()

    
    oco = adata.obs[column]
    if str(oco.dtype) in ['category', 'object']:
        # guess categorical
        uniq = len(oco.unique())
        if uniq > 15:
            dtype = 'skip'
        else:
            dtype = 'categorical'
    elif str(oco.dtype).startswith('int'):
        dtype = 'int'
    elif str(oco.dtype).startswith('float'):
        dtype = 'float'
    else:
        lg.warning(f"Unknown data type for {column} - {oco.dtype}")
        dtype = 'skip'

    nice = " ".join(
        [x.capitalize() for x in re.split(r'[\s_\_]', column)])
        
    return dict(name=column, type=dtype, nice=nice)


@meta.command('dump')
@click.argument('h5adfile')
def prepare_dump(h5adfile: str) -> None:
    "Dump metadata"
    import scanpy as sc
    
    lg.info("Loading data")
    adata = sc.read_h5ad(h5adfile, backed='r')

    if 'termite_obs_meta' in adata.uns:
        # see if the column is in the adata file
        print(adata.uns['termite_obs_meta'])
    else:
        print("No obs metadata in adata")


    
@meta.command('prepare')
@click.argument('h5adfile')
def prepare_meta(h5adfile: str) -> None:

    import scanpy as sc
    import pandas as pd

    basename = h5adfile.replace('.h5ad', '')

    lg.info("Loading data")
    adata = sc.read_h5ad(h5adfile)

    bom = None
    if 'termite_obs_meta' in adata.uns:
        # see if the column is in the adata file
        bom = adata.uns['termite_obs_meta']
        
    if not isinstance(bom, pd.DataFrame):
        bom = pd.DataFrame(columns = ['name', 'type', 'nice'])

    bom = bom.reset_index(drop=True)
    
    obs = adata.obs

    # start with DimRed columns
    lg.info("Processing obsm/dimred table")
    lg.info(f"  - keys: {adata.obsm_keys()}")
    for k in adata.obsm_keys():
        
        if k.startswith('_'):
            continue
        
        omat = adata.obsm[k]
        assert omat.shape[1] >= 2

        name0 = f"{k}/0"
        name1 = f"{k}/1"

        if name0 not in bom['name'].values:
            bom.loc[len(bom)] = dict(
                name=name0, type='dimred',
                nice=k.replace('X_', '').capitalize() + '/0')
        if name1 not in bom['name'].values:
            bom.loc[len(bom)] = dict(
                name=name1, type='dimred',
                nice=k.replace('X_', '').capitalize() + '/1')

    # regular obs columns
    lg.info("Processing obs table")
    for column in obs.keys():
        if column.startswith('_'):
            continue

        if column in bom['name'].values:
            continue

        bom.loc[len(bom)] = get_obs_column_metadata(adata, column)

    adata.uns['termite_obs_meta'] = bom
    lg.info("Saving files")
    adata.write(h5adfile)
    bom.to_csv(f"{basename}.obscol.tsv", sep="\t")
    

@meta.command('import')
@click.argument('h5adfile')
def import_meta(h5adfile: str) -> None:
    """Import metadata into the scanpy h5ad file."""
    
    import scanpy as sc
    import pandas as pd

    basename = h5adfile.replace('.h5ad', '')

    obs_meta_file = f'{basename}.obscol.tsv'
    if not os.path.exists(obs_meta_file):
        print(f"Did not find {obs_meta_file}")
        return
            
    obs_meta = pd.read_csv(obs_meta_file, sep="\t", index_col=0)
    
    adata = sc.read_h5ad(h5adfile) #/////backed='r+')
    adata.uns['termite_obs_meta'] = obs_meta
    adata.write(h5adfile)
    

@h5ad.command('import')
@click.argument('h5adfile')
@click.option('--expname')
@click.option('--datatype', default='raw', type=str)
@click.option('--layer', default=None, type=str)
def load_h5ad(h5adfile, expname, datatype, layer):

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
        elif dtype in ['int','float', 'dimred']:
            col['value'] = col['value'].astype(float)
            if db.table_exists(conn, 'obs_num'):
                conn.sql(f"""
                    DELETE FROM obs_num 
                     WHERE name='{colname}' 
                       AND experiment='{expname}'
                """)
            db.create_or_append(conn, 'obs_num', col)
        else:
            col['value'] = col['value'].astype(str)
            
            if db.table_exists(conn, 'obs_cat'):
                conn.sql(f"""
                    DELETE FROM obs_cat 
                     WHERE name='{colname}' 
                       AND experiment='{expname}'
                """)
            db.create_or_append(conn, 'obs_cat', col)


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
        dtype = get_obs_column_metadata(adata, colname)
        lg.info(f'processing obs | {dtype:10} | {colname}')

        total_obs_recs += len(col)
        _store_col(colname, col, dtype)

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
    if db.table_exists(conn, 'expr'):
        sql = f"""
            DELETE FROM expr 
             WHERE datatype='{datatype}' 
               AND experiment='{expname}'
            """
        conn.sql(sql)
        
    db.create_or_append(conn, 'expr', melted)
    
    db.status(conn)

