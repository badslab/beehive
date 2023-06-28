
from functools import lru_cache
import logging
import os
from pathlib import Path
import re
from typing import Dict

import click
import pandas as pd
import numpy as np

from termite import db, util, config


lg = logging.getLogger(__name__)


@click.group('h5ad')
def h5ad():
    pass


@lru_cache(16)
def get_metadata(h5adfile: str) -> Dict[str, pd.DataFrame]:

    # TODO: figure out why do I return the datafram in a dict?
    rv = {}
    
    obscol_mdfile = h5adfile.replace('.h5ad', '') + '.obscol.tsv'
    if os.path.exists(obscol_mdfile):
        rv['obscol'] = pd.read_csv(obscol_mdfile, sep="\t", index_col=0)
        
    return rv


@lru_cache(16)
def get_layerdata(h5adfile: str) -> pd.DataFrame:    
    layer_mdfile = h5adfile.replace('.h5ad', '') + '.datatypes.tsv'
    if os.path.exists(layer_mdfile):
        return pd.read_csv(layer_mdfile, 
                           sep="\t", index_col=0)
    else:
        raise FileNotFoundError(layer_mdfile)


    
COMMON_COLUMN_MAPPINGS = """
leiden_scVI | scVI cluster
""".strip().split("\n")
COMMON_COLUMN_MAPPINGS = \
    {a.strip():b.strip() for (a,b)
     in [x.split('|', 1) for x in COMMON_COLUMN_MAPPINGS]
     }

    
def get_obs_column_metadata(h5adfile, adata, column, force=False,
                            max_categories=30):
    """

    
    """

    all_md = get_metadata(h5adfile)
    md_obscol = all_md.get('obscol')

    # see if the column is in the obsocl metadata file
    if (not force) and (md_obscol is not None) and (column in md_obscol['name'].values):
        rv = md_obscol.loc[md_obscol['name'] == column].iloc[0].to_dict()
        return rv
        
    # nothing here? Guess the format from the raw data
    oco = adata.obs[column]
    example=",".join(map(str, oco.head(4)))
    
    if str(oco.dtype) in ['category', 'object', 'bool']:
        # guess categorical
        uniq = len(oco.unique())
        example = ",".join(map(str, oco.value_counts().sort_values()[:4].index))
        if uniq > max_categories:
            # too many categories
            dtype = 'skip'
        elif uniq == 1:
            # not very interesting either...
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

    if column in COMMON_COLUMN_MAPPINGS:
        nice = COMMON_COLUMN_MAPPINGS[column]
    else:
        nice = " ".join(
            [x.capitalize() for x in re.split(r'[\s_\-]', column)])
        
    return dict(name=column, type=dtype, nice=nice, example=example)        


@h5ad.group('prepare')
def prepare():
    pass


@prepare.command('experiment')
@click.option('-p', '--pubmed', default='unknown')
@click.option('-a', '--author', default='unknown')
@click.option('-t', '--title', default='unknown')
@click.option('-x', '--experiment')
@click.option('-y', '--year', type=int, default=1859)
@click.option('-f', '--force', type=bool, is_flag=True, default=False)
@click.option('-S', '--study')
@click.option('-v', '--version', type=int, default=1)
@click.option('-o', '--organism', default='unknown')
@click.argument('h5adfile')
def prepare_experiment(h5adfile: str,
                 author: str,
                 experiment: str,
                 force: bool,
                 year: int,
                 pubmed: str,
                 study: str,
                 version: int,
                 organism: str,
                 title: str) -> None:

    import scanpy as sc

    basename = h5adfile.replace('.h5ad', '')

    outfile = Path(f"{basename}.experiment.tsv")
    if outfile.exists() and not force:
        raise click.ClickException(
            f"{outfile} exists, use -f to overwrite")

    layerout = Path(f"{basename}.datatypes.tsv")
    if layerout.exists() and not force:
        raise click.ClickException(
            f"{layerout} exists, use -f to overwrite")

    
    lg.info(f"Loading data: {h5adfile}")
    adata = sc.read_h5ad(h5adfile)

    
    # find topgenes - to make suggestions in the interface
    c = adata.to_df()
    topgenes = c.sum()

    fgenes = topgenes.index.str.lower().isin(config.FAVORITE_GENES)
    topgenes.loc[fgenes] = topgenes.loc[fgenes] * 20

    nfgenes = topgenes.index.isin(config.NOT_FAVORITE_GENES)
    topgenes[nfgenes] = 0
    
    mtgenes = topgenes.index.str.lower().str.startswith('mt-')
    topgenes[mtgenes] = 0
    
    topgenes = topgenes[topgenes > 0]
    topgenes = topgenes.sort_values(ascending=False)
    topgenes = list(topgenes.head(10).index)
    topgenes=  ' '.join(topgenes)
    lg.info(f"Topgenes {topgenes}")
    

    # store data types
    datatypes = [
        dict(layer='X', datatype='Raw counts') ]
    
    for layer in adata.layers.keys():
        datatypes.append(
            dict(layer=layer, datatype=layer))
        
    datatypes = pd.DataFrame(datatypes)

    #
    # get dimred names
    #
    
    dimred_names = []
    for obsmname in adata.obsm_keys():
        if obsmname.startswith('_'):
            lg.warning(f"skipping obsm {obsmname}, starts with a '_'")
            continue
        
        obsm = adata.obsm[obsmname]

        if obsm.shape[1] < 2:
            lg.warning(f"skipping obsm {obsmname}, too few dimensions")
            #too few dimensions - skip
            continue

        if ' ' in obsmname:
            lg.warning(f"skipping obsm '{obsmname}', no spaces allowed")
            #too few dimensions - skip
            continue
        
        if not isinstance(obsm, np.ndarray):
            # this should be a numpy ndarray - if not skip
            lg.warning(f"skipping obsm {obsmname}, not a numpy.ndarray")
            continue

        dimred_names.append(obsmname)
    dimred_names = ' '.join(dimred_names)
    
    #
    # get experimental data
    #
    doi = 'unknown'
    pubmedfile = f"{basename}.pubmed.id"

    lg.info('pubmed ' + str(pubmed))
    if pubmed != 'unknown':
        lg.info(f'save pubmed {pubmed} for future reference')
        with open(pubmedfile, 'w') as F:
            F.write(pubmed)            
    elif os.path.exists(pubmedfile):
        pubmed = open(pubmedfile).read().strip()
        lg.info(f'reusing pubmed {pubmed} from file')

    if pubmed != 'unknown':            
        x = util.query_pubmed(pubmed)    
        if author == 'unknown':
            author = x['author']
        if title == 'unknown':
            title = x['title']
        if year == 1859:
            year = x['year']
        doi = x['doi']

    if experiment is None:
        experiment = os.path.basename(h5adfile).replace('.h5ad', '')
        
    if re.match(r'^[hmx]\.\w{3,10}\.\d{1,5}.*$', experiment):
        if organism == 'unknown':
            if experiment.startswith('h.'):
                organism = 'human'
            elif experiment.startswith('m.'):
                organism = 'mouse'
            elif experiment.startswith('x.'):
                organism = 'mixed'

        if study is None:
            study = experiment.split('.')[1]
            
    if study is None:
        study = experiment
        
    experiment_md = dict(
        experiment = experiment,
        author = author,
        study = study,
        version = version,
        title = title,
        doi = doi,
        year = year,
        pubmed = pubmed,
        organism = organism,
        dimred = dimred_names,
        topgenes = topgenes)

    experiment_md = pd.Series(experiment_md)
    experiment_md.to_csv(outfile, header=False, sep="\t")
    
    datatypes.to_csv(layerout, sep="\t")
    
    
@prepare.command('obs')
@click.option('-f', '--force', type=bool, is_flag=True, default=False)
@click.argument('h5adfile')
def prepare_obs(h5adfile: str, force: bool) -> None:

    import scanpy as sc

    basename = h5adfile.replace('.h5ad', '')

    outfile = Path(f"{basename}.obscol.tsv")
    if outfile.exists() and not force:
        raise click.ClickException(
            f"{outfile} exists, use -f to overwrite")
        
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
            = get_obs_column_metadata(h5adfile, adata, column, force=force)

    lg.info("Saving MD files")
    md_obscol.to_csv(outfile, sep="\t")


@h5ad.command('import')
@click.argument('h5adfile')
@click.option('-f', '--forget', type=bool, is_flag=True, default=False)
@click.option('--expname')
@click.option('--datatype', default='raw', type=str)
@click.option('--layer', default=None, type=str)
def load_import(h5adfile, expname, datatype, layer, forget):

    import scanpy as sc
    import pandas as pd
    conn = db.get_conn()

    basename = h5adfile.replace('.h5ad', '')
    
    # start by reading the experimental metadata
    expdata = pd.read_csv(basename + '.experiment.tsv', sep="\t",
                          index_col=0, header=None).T

    expname = expdata.loc[1]['experiment']

    if forget:
        lg.info(f"Forgetting about {expname}")
        db.forget(expname)
        
    lg.info(f"import experiment {expname}")

    if pd.isna(expdata.loc[1]['dimred']):
        expdata.loc[1]['dimred'] = ''
        
    db.create_or_append('experiment_md', expdata)

    adata = sc.read_h5ad(h5adfile)
        
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
        else:  # assuming categorical
            from pandas.api.types import is_float_dtype
            if is_float_dtype(col['value']):
                try:
                    # try to see if these really are integers..
                    print('attempt cast to integer first')
                    newcol = col['value'].astype(int)
                    if abs(newcol - col['value']).sum() == 0:
                        # it is an int:
                        print('success!')
                        col['value'] == newcol
                except:
                    print('failed conversion to int, keep float')
                
                
            col['value'] = col['value'].astype(str)

            
            if db.table_exists('obs_cat'):
                conn.sql(f"""
                    DELETE FROM obs_cat 
                     WHERE name='{colname}' 
                       AND experiment='{expname}'
                """)
            db.create_or_append('obs_cat', col)


    dimred_names = expdata.loc[1].fillna('')['dimred'].split()
    for obsmname in dimred_names:

        obsm = adata.obsm[obsmname]
        lg.info(f"storing DimRed {obsmname} {obsm.shape} {type(obsm)}")
        obsmdata = pd.DataFrame(obsm[:, :2]).astype(float)
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

    layerdata = get_layerdata(h5adfile)
    for _, layer in layerdata.iterrows():
        lg.info(f"Processing layer {layer['layer']} / {layer['datatype']}")
        if layer['datatype'] == 'skip':
            continue

        if layer['layer'] == 'X':
            x = adata.to_df()
        else:
            x = adata.to_df(layer=layer['layer'])

        # chunk this - using too much memory:
        
        x.index.name = 'obs'
        x.columns.name = 'gene'

        
        #remove old data
        lg.info("remove old data")
        if db.table_exists('expr'):
            sql = f"""
                DELETE FROM expr 
                 WHERE datatype='{datatype}' 
                   AND experiment='{expname}'
            """
            conn.sql(sql)

        chunksize = 5000
        for ic in range(0, x.shape[0], chunksize):
            chunk = x.iloc[ic:ic+chunksize,:]
            melted = chunk.reset_index().melt(id_vars='obs')
            melted['experiment'] = expname
            melted['datatype'] = layer['datatype']
            lg.info(f"chunk {ic}/{x.shape[0]} - melt {melted.shape[0]:_d}")
            db.create_or_append('expr', melted)

    tablecount = db.all_table_count()
    for t in sorted(tablecount):
        c = tablecount[t]
        print(f"{t:<20s} : {c:>14_d}")

