
import logging
from functools import partial

import click
from click.core import Context

from . import db
from . import metadata_tools as mdtools
from . import util

lg = logging.getLogger(__name__)


@click.command("upload")
@click.argument('h5ad', type=click.Path(exists=True))
@click.option("-c", "skip_counts", is_flag=True, default=False,
              help="Skip count table.")
@click.option("-o", "skip_obs", is_flag=True, default=False,
              help="Skip obs table.")
@click.option("-m", "skip_obsm", is_flag=True, default=False,
              help="Skip obsm table.")
@click.pass_context
def upload(ctx: Context,
           h5ad: str,
           skip_counts: bool,
           skip_obs: bool,
           skip_obsm: bool,) -> None:
    """Upload an h5ad file to the database."""

    # late import to speed up matters
    import pandas as pd
    import scanpy as sc

    # database object.
    chdb = ctx.obj['chdb']

    # import the h5ad file
    adata = sc.read_h5ad(h5ad)

    # check if there are problems.
    problems = mdtools.check_2(adata)
    if problems:
        for p in problems:
            print(p)
        return

    # helper fuction
    mdget = partial(util.mdget, data=adata.uns['cellhive'])

    # preparing the basis for the expdata table
    # it's however quite denormalized - one row per
    # layer.
    expdata = mdget('metadata')

    study = mdget('metadata', 'study')
    version = mdget('metadata', 'version', default='0')
    expname = mdget('metadata', 'experiment')


    expdata['study_id'] = chdb.get_id('experiment_md', 'study', study)

    expdata['full_experiment'] = full_exp_name = \
        f"{study}__{expname}__{version}"

    expdata['full_experiment_id'] = \
        chdb.get_id('experiment_md', 'full_experiment', full_exp_name)

    #OBSM
    if not skip_obsm:
        obsm_data = mdtools.obsm(adata)
        assert isinstance(obsm_data, pd.DataFrame)
        for obsm_name in obsm_data:
            obsm = obsm_data[obsm_name]
            if str(obsm['ignore']) == 'True':
                continue
            lg.info(f'import obsm {obsm_name}')
            obd = adata.obsm[str(obsm_name)]
            lg.info(f"storing obsm col {obsm_name}")
            chdb.store_obscol(
                col=pd.Series(obd[:,0], index=adata.obs_names),
                name=f"{obsm_name}/0",
                dtype='float',
                exp_id=expdata['full_experiment_id'])

            chdb.store_obscol(
                col=pd.Series(obd[:,1], index=adata.obs_names),
                name=f"{obsm_name}/1",
                dtype='float',
                exp_id=expdata['full_experiment_id'])

    #OBS
    if not skip_obs:
        obs_data = mdtools.obs(adata)
        for _, obs in obs_data.iterrows():
            od = obs.to_dict()
            if od.get('ignore') is True:
                lg.info(f"Ignoring obs col {obs['name']}")
                continue

            lg.info(f"storing obs col {obs['name']}")
            chdb.store_obscol(
                col=adata.obs[od['name']],
                name=obs['name'],
                dtype=obs['dtype'],
                exp_id=expdata['full_experiment_id'])

    # Layers!
    layerdata = mdtools.layers(adata)
    assert layerdata is not None
    for _, linfo in layerdata.items():
        if str(linfo['ignore']) == 'True':
            continue

        expdata_l = expdata.copy()
        expdata_l['layer_name'] = lname = linfo['name']
        expdata_l['layer_type']  = linfo['layer_type']

        dataset = f"{full_exp_name}__{lname}"
        expdata_l['dataset'] = dataset
        expdata_l['dataset_id'] = dataset_id = \
            chdb.get_id('experiment_md', 'dataset', dataset)

        lg.info("Store layer info")
        chdb.uac_experiment_md(expdata_l)


        if not skip_counts:
            lg.info("Start count import")

            chdb.import_count_table(
                dataset_id=dataset_id,
                adata=adata,
                layer=lname)
