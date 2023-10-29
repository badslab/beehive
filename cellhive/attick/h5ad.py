

import logging
import os
import re
from array import array
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional

import click
import duckdb
import numpy as np
import pandas as pd
import scanpy as sc
from rich.pretty import pprint
from scanpy import AnnData

lg = logging.getLogger(__name__)

# helper utilities to get/set metadata in a scanpy h5ad file



def remove(adata: AnnData,
           key: str):
    if 'termite' not in adata.uns_keys():
        return
    if key not in adata.uns['termite']['metadata']:
        return
    del adata.uns['termite']['metadata'][key]


def set1(adata: AnnData,
         key: str,
         value: Any) -> None:
    check(adata)
    adata.uns['termite']['metadata'][key] = value


COMMON_COLUMN_MAPPINGS = """
leiden_scVI | scVI cluster
""".strip().split("\n")
COMMON_COLUMN_MAPPINGS = \
    {a.strip():b.strip() for (a,b)
     in [x.split('|', 1) for x in COMMON_COLUMN_MAPPINGS]
     }


def get_obs_column_metadata(adata, column,
                            max_categories=30):
    """
    # determine what datatye a column is
    """

    # Guess the format from the raw data
    oco = adata.obs[column]
    no_uniq = len(oco.unique())
    example=",".join(map(str, oco.head(4)))

    forced = False
    tdata = adata.uns['termite']
    #if column == 'batch':
    #    print(tdata['obs_force'])
    if column in  tdata['obs_force']:
        odtype = tdata['obs_force'][column]
        forced = True
    else:
        odtype = str(oco.dtype)

    if odtype in ['categorical', 'category', 'object', 'bool']:
        # guess categorical
        uniq = len(oco.unique())
        example = ",".join(map(str, oco.value_counts().sort_values()[:4].index))
        dtype = 'categorical'
        # unless....
        if not forced:
            if uniq > max_categories:
                # too many categories
                dtype = 'skip'
            elif uniq == 1:
                # not very interesting either...
                dtype = 'skip'
    elif 'int' in odtype:
        dtype = 'int'
    elif odtype.startswith('float'):
        dtype = 'float'
    elif odtype == 'skip':
        dtype = 'skip'
    else:
        lg.warning(f"Unknown data type for {column} - {oco.dtype}")
        dtype = 'skip'

    if column in tdata['obs_nice']:
        nice = tdata['obs_nice'][column]
    elif column in COMMON_COLUMN_MAPPINGS:
        nice = COMMON_COLUMN_MAPPINGS[column]
    else:
        nice = (" ".join(
            [x.capitalize() for x in re.split(r'[\s_\-\.]+', column)])
                ).strip()
    return dict(name=column, dtype=dtype,
                nice_name=nice, no_uniq=no_uniq,
                example=example)



# @click.argument('h5adfile')
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
    datatypes = ['X']
    for layer in adata.layers.keys():
        datatypes.append(layer)

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
        h5ad_layers = ' '.join(datatypes),
        layer = datatypes[0],
        datatype = 'counts',
        organism = organism,
        dimred = dimred_names,
        topgenes = topgenes)

    experiment_md = pd.Series(experiment_md)
    experiment_md.to_csv(outfile, header=False, sep="\t")


def prepare_obs(adata: AnnData) -> pd.DataFrame:

    md_obscol = {}
    obs = adata.obs

    # regular obs columns
    lg.info("Processing obs table")

    for column in obs.keys():
        if column.startswith('_'):
            continue
        md_obscol[column] \
            = get_obs_column_metadata(adata, column)

    return md_obscol


def preprocess_catcol(col):
    """ Prepare a column called as categorical

    is it float? try to cast to int first
    convert to string

    """

    from pandas.api.types import is_float_dtype
    if is_float_dtype(col):
        try:
            # try to see if these really are integers..
            lg.debug('Attempt cast to integer first')
            newcol = col.astype(int)
            if abs(newcol - col).sum() == 0:
                # it is an int:
                lg.debug('success casting to int')
                return newcol
        except:
            print('failed conversion to int, keep float')
    return col.astype(str)


def import_experiment_md(adata: AnnData,
                         layername: str,
                         layertype: str,
                         chunksize: int = 5000) -> int:

    expdata = deepcopy(adata.uns['termite']['metadata'])

    expdata['layername'] = layername
    expdata['layertype'] = layertype
    experiment = expdata['experiment']
    dataset = f"{experiment}.{layername}"
    expdata['dataset'] = dataset

    if 'dimred' not in expdata:
        expdata['dimred'] = ''

    exp_id = db.autoincrementor('dataset_md', 'experiment', experiment)
    dataset_id = db.autoincrementor('dataset_md', 'dataset', dataset)

    lg.info(f"using experiment {experiment} {exp_id}")
    lg.info(f"      dataset    {dataset} {dataset_id}")

    expdata['experiment_id'] = exp_id
    expdata['dataset_id'] = dataset_id

    # remove if there is already a dataset record
    try:
        db.raw_sql(f'DELETE FROM dataset_md WHERE dataset_id={dataset_id}')
    except duckdb.CatalogException:
        pass  # table does not exist?

    for f in 'abstract author doi title organism'.split():
        if f not in expdata:
            expdata[f] = ''
        else:
            expdata[f] = str(expdata[f])

    for f in 'pubmed year'.split():
        if f not in expdata:
            expdata[f] = 0
        else:
            expdata[f] = int(expdata[f])

    expdata = pd.DataFrame([expdata])

    db.create_or_append('dataset_md', expdata)

    return dataset_id


def store_one_obs_col(exp_id: int,
                      colname: str,
                      col: pd.DataFrame,
                      dtype: str,
                      original_name: str,
                      dimred_name: str='',
                      dimred_dim: int=-1):

    conn = db.get_conn()
    col.columns = ['value']

    col.index.name = 'obs'
    col = col.reset_index()
    col['name'] = colname
    col['exp_id'] = exp_id

    colmeta = dict(
        name=colname,
        exp_id=exp_id,
        dtype=dtype,
        original_name=original_name,
        dimred_dim=dimred_dim,
        dimred_name=dimred_name,
        no_cat=-1,
        )

    if dtype == 'skip':
        pass
    elif dtype in ['int', 'float', 'dimred']:
        col['value'] = col['value'].astype(float)
        if db.table_exists('obs_num'):
            conn.sql(f"""
                DELETE FROM obs_num
                 WHERE name='{colname}'
                   AND exp_id='{exp_id}' """)
        db.create_or_append('obs_num', col)

    else:  # assuming categorical
        col['value'] = preprocess_catcol(col['value'])
        colmeta['no_cat'] = len(col['value'].unique())

        if db.table_exists('obs_cat'):
            conn.sql(f"""
                DELETE FROM obs_cat
                 WHERE name='{colname}'
                   AND exp_id='{exp_id}'
            """)
        db.create_or_append('obs_cat', col)

    try:
        db.raw_sql(
            f""" DELETE FROM help_obs
                  WHERE exp_id={exp_id}
                    AND name='{colname}' """)
    except duckdb.CatalogException:
        pass # db does not exist...

    db.create_or_append('help_obs', pd.Series(colmeta))


def import_experiment_obs(exp_id: int,
                          adata: AnnData):

    lg.info(f"exp id    : {exp_id}")

    total_obs_recs = 0

    for colname in adata.obs:
        if colname.startswith('_'):
            continue
        col = adata.obs[[colname]].copy()
        colinfo = get_obs_column_metadata(adata, colname)
        if colinfo['dtype'] == 'skip':
            continue

        nice = str(colinfo['nice_name'])
        dtype = str(colinfo['dtype'])
        lg.info(f'processing obs | {dtype:14s} | {colname}')
        total_obs_recs += len(col)
        store_one_obs_col(
            exp_id=exp_id, colname=nice, col=col,
            dtype=dtype, original_name=colname)

    lg.info(f"total obs records {total_obs_recs}")


def import_experiment_dimred(exp_id: int,
                             adata: AnnData):

    lg.info(f"exp id    : {exp_id}")

    total_obs_recs = 0

    dimred_names = expdata.fillna('')['dimred'].split()
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


def import_counts(dataset_id: int,
                  adata: AnnData,
                  layer: str,
                  normalize: Optional[str] = None,
                  chunksize: int = 50000,
                  ) -> None:

    conn = db.get_conn()
    lg.info("Start storing expression matrix")
    lg.info(f"Processing layer {layer}")

    if layer == 'X':
        x = adata.to_df()
    else:
        x = adata.to_df(layer=layer)

    # chunk this - using too much memory:
    x.index.name = 'obs'
    x.columns.name = 'gene'

    if normalize == 'logrpm':
        x = np.log1p(10000 * x.copy().divide(x.sum(1), axis=0))


    #remove old data
    lg.info("remove old data")
    if db.table_exists('expr'):
        sql = f"""
            DELETE FROM expr
             WHERE dataset_id={dataset_id}"""
        conn.sql(sql)


    lg.info("start expression data upload")
    for ic in range(0, x.shape[0], chunksize):
        chunk = x.iloc[ic:ic+chunksize,:]
        # chunk_rnk = rnk.iloc[ic:ic+chunksize,:]

        melted = chunk.reset_index().melt(id_vars='obs')
        melted['value'] = melted['value'].astype(float)  # ensure!
        melted['dataset_id'] = dataset_id

        # melted_rnk = chunk.reset_index().melt(id_vars='obs', value_name='rank')
        # melted_rnk['rank'] = melted_rnk['rank'].astype(float)  # ensure!

        # to be sure!
        # assert melted[['obs', 'gene']].equals(melted_rnk[['obs', 'gene']])
        # melted['rank'] = melted_rnk['rank']

        lg.info(f"chunk {ic}/{x.shape[0]} - melt {melted.shape[0]:_d}")
        db.create_or_append('expr', melted)

    # ensure index
    # print(db.raw_sql('create index idx_expr_eg on expr (exp_id, gene)'))


#@click.command('import')
#@click.argument('h5adfile')
#@click.option('-f', '--forget', type=bool, is_flag=True, default=False)
#@click.option('--expname')
#@click.option('--datatype', default='raw', type=str)
#@click.option('--chunksize', default=20000, type=int)
#@click.option('--layer', default=None, type=str)
def h5ad_import_old(h5adfile, expname, datatype, layer, chunksize, forget):

    import pandas as pd
    import scanpy as sc
    conn = db.get_conn()

    basename = h5adfile.replace('.h5ad', '')

    # start by reading the experimental metadata
    expdata = pd.read_csv(basename + '.experiment.tsv', sep="\t",
                          index_col=0, header=None).T.loc[1]

    expname = expdata.loc['experiment']

    if forget:
        lg.info(f"Forgetting about {expname}")
        db.forget(expname)

    lg.info(f"import experiment {expname}")

    if pd.isna(expdata.loc['dimred']):
        expdata.loc['dimred'] = ''

    exp_id = db.get_experiment_id(expname)


    if exp_id is not None:
        # a previous record exists - removing
        sql = f"""DELETE FROM experiment_md
                   WHERE exp_id = {exp_id}"""
        db.raw_sql(sql)
    else:
        # no prev. record - caluclate a good new exp_id
        try:
            max_exp = db.raw_sql(
                '''SELECT MAX(exp_id) AS exp_id
                     FROM experiment_md''').iloc[0,0]
        except duckdb.CatalogException:
            # assuming the table does not exist yet...
            max_exp = 0
        exp_id = max_exp + 1

    expdata['exp_id'] = exp_id

    db.create_or_append('experiment_md', expdata)

    adata = sc.read_h5ad(h5adfile)

    #print(adata.var.head(3).T)
    #exit()

    # remove obs_names - enforce integers to reduce database size
    adata.obs_names = list(range(len(adata)))

    lg.info(f"data name : {expname}")
    lg.info(f"exp id    : {exp_id}")
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
        col['exp_id'] = exp_id

        if dtype == 'skip':
            pass
        elif dtype in ['int', 'float', 'dimred']:

            col['value'] = col['value'].astype(float)
            if db.table_exists('obs_num'):
                conn.sql(f"""
                    DELETE FROM obs_num
                     WHERE name='{colname}'
                       AND exp_id='{exp_id}' """)
            db.create_or_append('obs_num', col)

        else:  # assuming categorical
            col['value'] = preprocess_catcol(col['value'])

            if db.table_exists('obs_cat'):
                conn.sql(f"""
                    DELETE FROM obs_cat
                     WHERE name='{colname}'
                       AND exp_id='{exp_id}'
                """)
            db.create_or_append('obs_cat', col)

    dimred_names = expdata.fillna('')['dimred'].split()
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

    layer = expdata['layer']
    lg.info(f"Processing layer {layer}")

    if layer == 'X':
        x = adata.to_df()
    else:
        x = adata.to_df(layer=layer)


    # chunk this - using too much memory:

    x.index.name = 'obs'
    x.columns.name = 'gene'

    lg.info("calculate ranks")
    rnk = x.rank(pct=True, method='min', na_option='bottom')

    #remove old data
    lg.info("remove old data")
    if db.table_exists('expr'):
        sql = f"""
            DELETE FROM expr
             WHERE exp_id={exp_id}"""
        conn.sql(sql)


    lg.info("start expression data upload")
    for ic in range(0, x.shape[0], chunksize):
        chunk = x.iloc[ic:ic+chunksize,:]
        chunk_rnk = rnk.iloc[ic:ic+chunksize,:]

        melted = chunk.reset_index().melt(id_vars='obs')
        melted['value'] = melted['value'].astype(float)  # ensure!
        melted['exp_id'] = exp_id

        melted_rnk = chunk.reset_index().melt(id_vars='obs', value_name='rank')
        melted_rnk['rank'] = melted_rnk['rank'].astype(float)  # ensure!

        # to be sure!
        assert melted[['obs', 'gene']].equals(melted_rnk[['obs', 'gene']])
        melted['rank'] = melted_rnk['rank']

        lg.info(f"chunk {ic}/{x.shape[0]} - melt {melted.shape[0]:_d}")
        db.create_or_append('expr', melted)

    # ensure index
    # print(db.raw_sql('create index idx_expr_eg on expr (exp_id, gene)'))

