
import logging
import os
from pathlib import Path
from typing import Optional
import textwrap

import anndata
import click
import pandas as pd
from rich.logging import RichHandler
import scanpy as sc

import termite.cli.db
import termite.h5ad
import termite.db
import termite.util
from termite.globals import globals


FORMAT = "%(message)s"
logging.basicConfig(
    level=logging.INFO, format=FORMAT, datefmt="[%X]",
    handlers=[RichHandler()])


lg = logging.getLogger('termite')
lg.info("Start termite!")


def in_jupyter() -> bool:
    """Check if we're called from jupyter."""
    return 'JPY_PARENT_PID' in os.environ


def where_is_adata(func):
    "Decorator to get adata from *args or `globals`"
    def call(*args, **kwargs):
        if len(args) == 1 and isinstance(args[0], sc.AnnData):
            adata = args[0]
        else:
            adata = globals.adata
        return func(adata, **kwargs)
    return call


@click.group("cli")
def cli() -> None:
    pass


@where_is_adata
def md(adata, **kwargs):
    """Show or set metadata"""

    if 'pubmed' in kwargs:
        pubmed_id = kwargs['pubmed']
        pmd = termite.util.query_pubmed(pubmed_id)
        termite.h5ad.set1(adata, 'author', pmd['author'])
        termite.h5ad.set1(adata, 'title', pmd['title'])
        termite.h5ad.set1(adata, 'year', pmd['year'])
        termite.h5ad.set1(adata, 'abstract', pmd['abstract'])
        termite.h5ad.set1(adata, 'doi', pmd['doi'])

    for k, v in kwargs.items():
        termite.h5ad.set1(adata, k, v)

    if 'termite' not in adata.uns:
        print("no annotations")
        return

    if in_jupyter():
        return adata.uns['termite']['metadata']

    # not in jupyter - print to screen
    for k, v in adata.uns['termite']['metadata'].items():
        print(k + ": " + "\n".join(
            textwrap.wrap(str(v), subsequent_indent='...')))


@where_is_adata
def get_layerdata(adata, *args):
    """Create stats on the layers in this adata"""
    ldata = adata.uns['termite']['layers']
    d = []
    def check_layer(name: str,
                    X):

        if isinstance(
                X, anndata._core.sparse_dataset.SparseDataset):
            X = X.value.todense()

        row = {}
        row['name'] = name
        row['dtype'] = X.dtype
        row['rows'] = X.shape[0]
        row['columns'] = X.shape[1]
        row['entries'] = entries = X.shape[0] * X.shape[1]
        row['min'] = X.min()
        row['max'] = X.max()

        zeros = len((X==0).nonzero()[0])
        perc_zeros = 100 * zeros / entries
        row['zeros'] = zeros
        row['% zeros'] = perc_zeros
        d.append(row)
        row['termite_type'] = ldata.get(name, {}).get('type', '-')
        row['import?'] = ldata.get(name, {}).get('load', '-')

    check_layer('X', adata.X)
    for layer in adata.layers.keys():
        check_layer(layer, adata.layers[layer])

    df = pd.DataFrame(d)
    return df


@where_is_adata
def obsm(adata: sc.AnnData,
         name: Optional[str] = None,
         load: bool = True):

    termite.h5ad.check(adata)
    ao = adata.uns['termite']['obsm']

    if name is not None:
        assert name in adata.obsm_keys()
        if name not in ao:
            ao[name] = {}
        if load:
            ao[name]['load'] = True

    rv = []
    for o in adata.obsm_keys():
        rv.append(dict(name=o,
                       dim=adata.obsm[o].shape[1],
                       load=ao.get(o, {}).get('load', False),
                       ))
    rv = pd.DataFrame(rv).set_index('name').T
    return rv


@where_is_adata
def layers(adata, name=None, ltype=None, load=None) -> Optional[pd.DataFrame]:
    """Get or set layer data.

    Usage:
       layers(): returns a dataframe
       layers(name='layername', [type='layertype', [load=True/False]]):
    """
    if name is not None:
        ldata = adata.uns['termite']['layers']
        if name not in ldata:
            ldata[name] = {}

        if ltype is not None:
            ldata[name]['type'] = ltype
        if load is not None:
            ldata[name]['load'] = load
        return

    df = get_layerdata()
    return df.T


@where_is_adata
def save(adata,
         filename: Optional[str] = None):

    if filename is None:
        h5adfile = h5adfile
    else:
        h5adfile = Path(filename)

    lg.info(f"Saving adata to {h5adfile}")
    adata.write_h5ad(h5adfile)


@where_is_adata
def check(adata: sc.AnnData):
    if not globals.init:
        lg.warning("no adata loaded?")
    else:
        termite.h5ad.check(adata)


def register(h5adfile: Path | str,
             adata) -> None:

    if isinstance(h5adfile, str):
        h5adfile = Path(h5adfile)

    basename = h5adfile.name
    if not basename.endswith('.h5ad'):
        lg.warning("h5adfile must end with '.h5ad' extension")
    experiment = basename[:-5]
    lg.info(f"experiment name: {experiment}")

    globals.h5adfile = h5adfile
    globals.adata = adata
    globals.init = True
    check()

    termite.h5ad.set1(globals.adata, 'experiment', experiment)
    if not 'study' in globals.adata.uns['termite']['metadata']:
        globals.adata.uns['termite']['metadata']['study'] = experiment



def load(filename: Optional[str] = None) -> None:
    if filename is None:
        if not globals.init:
            lg.warning("No h5adfile loaded")
            return
        h5adfile = globals.h5adfile
    else:
        h5adfile = Path(filename)

    lg.info(f"Load h5ad file {h5adfile} into adata")
    basename = h5adfile.name
    if not basename.endswith('.h5ad'):
        lg.warning("h5adfile must end with '.h5ad' extension")
    experiment = basename[:-5]
    lg.info(f"experiment name: {experiment}")
    globals.h5adfile = h5adfile
    globals.adata = sc.read_h5ad(h5adfile)
    globals.init = True
    check()
    termite.h5ad.set1(globals.adata, 'experiment', experiment)
    if 'study' not in globals.adata.uns['termite']['metadata']:
        globals.adata.uns['termite']['metadata']['study'] = experiment


def obs_conv_int(name):
    globals.adata.obs[name] = globals.adata.obs[name].astype(int)
def obs_conv_float(name):
    globals.adata.obs[name] = globals.adata.obs[name].astype(float)
def obs_conv_skip(name):
    pass
def obs_conv_categorical(name):
    globals.adata.obs[name] = globals.adata.obs[name].astype(str)


OBS_CONV_FUNCTIONS = {
    'int': obs_conv_int,
    'float': obs_conv_float,
    'categorical': obs_conv_categorical,
    'skip': obs_conv_skip}


def obs(name=None, dtype=None, nice_name=None, what=None) -> Optional[pd.DataFrame]:
    termite.h5ad.check(globals.adata)
    if name is not None:
        td = globals.adata.uns['termite']
        if dtype is not None:
            if dtype not in ['int', 'float', 'categorical', 'skip']:
                lg.error("dtype must be one of: int, float, categorical or skip")
                return
            else:
                lg.warning(f"Convert obs column {name} to {dtype}")
                OBS_CONV_FUNCTIONS[dtype](name)
                td['obs_force'][name] = dtype

        if nice_name is not None:
            td['obs_nice'][name] = nice_name

    obs_data = pd.DataFrame(termite.h5ad.prepare_obs(globals.adata)).T
    obs_data = obs_data.reset_index(drop=True)
    if name is not None:
        return obs_data[obs_data['name'] == name].T
    else:
        if what == 'all':
            pass
        elif what is not None:
            obs_data = obs_data[obs_data['dtype'] == what]
        else:
            obs_data = obs_data[obs_data['dtype'] != "skip"]
        return obs_data


def todb(logrpm: bool = True,
         skip_layers: bool = False,
         skip_obs: bool = False,
         skip_obsm: bool = False,
         dim2load: int=2 ) -> None:

    adata = globals.adata
    termite.h5ad.check(adata)

    # remove obs_names - enforce numbers to reduce database size
    expname = adata.uns['termite']['metadata']['experiment']

    exp_id = termite.db.autoincrementor('dataset_md', 'experiment', expname)
    lg.info(f"Storing experiment {expname} with id {exp_id}")

    if not skip_obs:
        termite.h5ad.import_experiment_obs(exp_id, adata)

    if not skip_obsm:
        for (o, od) in adata.uns['termite']['obsm'].items():
            if not od.get('load'):
                continue
            lg.info(f"loading obsm: {o}")
            # default - load first 2 dimensions
            obsm = adata.obsm[o]
            for i in range(min(dim2load, obsm.shape[1])):
                c = pd.DataFrame(
                    pd.Series(adata.obsm[o][:,i],
                              index=adata.obs_names))
                colname = f"{o}/{i:02d}"
                termite.h5ad.store_one_obs_col(
                    exp_id=exp_id,
                    colname=colname,
                    col=c,
                    original_name='-',
                    dtype='dimred',
                    dimred_name=o,
                    dimred_dim=i)


    layerdata = adata.uns['termite']['layers']

    logrpm_in_adata = False
    for lname, ldata in layerdata.items():
        if ldata.get('type') == 'logrpm':
            logrpm_in_adata = True

    for lname, ldata in layerdata.items():
        if not ldata.get('load'):
            lg.info(f"Skipping load of layer {lname}")
        ltype = ldata['type']
        if not ltype:
            lg.error(f"cannot load layer {lname} "
                     + "with no type annotated")
            exit(-1)

        dataset_id = termite.h5ad.import_experiment_md(
            adata, lname, ltype)

        if not skip_layers:
            termite.h5ad.import_counts(
                dataset_id, adata, lname)


        if ltype == 'raw' and (not logrpm_in_adata) and logrpm:
            dataset_id_2 = termite.h5ad.import_experiment_md(
                adata, 'auto_logrpm', 'logrpm')
            if not skip_layers:
                termite.h5ad.import_counts(
                    dataset_id_2, adata, lname, normalize='logrpm')


def test():
    # CREATE TABLE help_gene AS
    sql = """
         select hg.dataset_id, hg.gene, hg.sumval,
                (select hg.sumval / quantile_cont(sumval, 0.99) from help_gene
                  where dataset_id=hg.dataset_id) as relsumval
           from help_gene as hg
          limit 20
          """
    return termite.db.raw_sql(sql)



def sql(sql):
    return termite.db.raw_sql(sql)

def helptables():
    lg.warning("creating support tables")
    termite.db.helptables()

@cli.command("todb")
@click.option("-l", "--logrpm", is_flag=True,
              default=True,
              help="autoconvert the raw count layer to logrpm")
@click.option("-o", "--skipobs", is_flag=True,
              default=False, help="Skip loading obs data")
@click.option("-s", "--skiplayers", is_flag=True,
              default=False, help="Skip loading all count data")
@click.argument('h5adfile')
def cli_todb(h5adfile, logrpm, skipobs, skiplayers):
    load(h5adfile)
    todb(logrpm=logrpm, skip_layers=skiplayers, skip_obs=skipobs)


@cli.command()
@click.argument("h5adfile")
def shell(h5adfile):
    load(h5adfile)
    lg.info("Dropping into ipython")
    from IPython import embed
    embed()



#cli.add_command(termite.h5ad.prepare)
#cli.add_command(termite.h5ad.h5ad_import)
#cli.add_command(termite.diffexp.de_run)
cli.add_command(termite.cli.db.db_group)


def run():
    cli()
