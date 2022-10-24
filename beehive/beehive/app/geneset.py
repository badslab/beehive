"""helpers for the gene expression app."""

from functools import partial
import copy
import hashlib
import logging
import os
import pickle
import sqlite3
import time

import pandas as pd
from pathlib import Path
from pprint import pprint
import gseapy as gp
import polars as pl
import shutil
import yaml

import typer
from typer import echo
from typing import List, Optional, Dict

import beehive
from beehive import util, expset
from beehive.util import dict_set, diskcache, query_pubmed


app = typer.Typer()

lg = logging.getLogger(__name__)


def load_yaml(filename: Path):
    with open(filename) as F:
        d = yaml.load(F, Loader=yaml.SafeLoader)
    return d


def get_hash(*args, **kwargs):
    h = hashlib.sha512()
    for a in sorted([str(x).lower() for x in args]):
        h.update(a.encode())

    for k, v in sorted(kwargs.items()):
        if isinstance(v, (tuple, list, dict)):
            continue
        h.update(str(k).lower().encode())
        h.update(str(v).lower().encode())
    return h.hexdigest()


def get_geneset_db():
    geneset_db_folder = beehive.BASEDIR / 'geneset'

    if not geneset_db_folder.exists():
        geneset_db_folder.mkdir(parents=True)

    return sqlite3.connect(
        geneset_db_folder / 'geneset_db.sqlite')


def get_geneset_groups():
    gsdb = get_geneset_db()
    groups = pd.read_sql('SELECT * FROM groups', gsdb)
    yield from groups.iterrows()


def get_geneset_genes(group_hash):
    gsdb = get_geneset_db()
    gsets = pd.read_sql(
        f'''SELECT * FROM genesets
            WHERE group_hash = '{group_hash}'
        ''', gsdb)
    return gsets


def run_one_gsea(args):

    group_hash, lfc_col, rnk, genedict = args

    results = gp.prerank(
        rnk=rnk, gene_sets=genedict,
        threads=3, min_size=5,
        max_size=1000, permutation_num=5000,
        outdir=None, seed=6,
        verbose=True)

    res2 = results.res2d[['Term', 'NES', 'FDR q-val']]
    res2.columns = ['set_hash', 'nes', 'fdr']
    res2.set_index('set_hash')

    lg.info(f"finished one gsea {group_hash} - {lfc_col}")
    return group_hash, lfc_col, res2


def run_one_gsea_cached(args):
    uid = util.UID(*args, length=12)
    cache_folder = beehive.BASEDIR / 'gsea' / 'cache'
    
    if not cache_folder.exists():
        try:
            cache_folder.mkdir(parents=True)
        except FileExistsError:
            pass

    cache_file = cache_folder / uid

    if cache_file.exists():
        with open(cache_file, 'rb') as F:
            lg.info(f"loading from cache: {uid}")
            rv = pickle.load(F)
    else:
        rv = run_one_gsea(args)
        lg.info(f"saving to cache: {uid}")
        with open(cache_file, 'wb') as F:
            pickle.dump(rv, F)

    return rv


def get_geneset_dicts() -> dict:
    """Return a dict of dictionaries
    """

    rv = {}

    for _, g in get_geneset_groups():
        gsets = get_geneset_genes(g['group_hash'])
        lg.debug(
            f"    | {g['study_title'][:40]} "
            f"| {g['group_title'][:40]}")

        def gg(x):
            return list(set(x.split()))

        genedict = {
            row['geneset_hash']: list(sorted(gg(row['genes'])))
            for (_, row) in gsets.iterrows()}

        rv[g['group_hash']] = genedict

    return rv


@ app.command('gsea')
def gsea(
        dsid: str = typer.Argument(..., help='Dataset'),
        threads: int = typer.Option(32, '-j', '--threads'), ):

    from multiprocessing import Pool

    output_file = util.find_prq(dsid, 'gsea', check_exists=False)
    lg.info(f'Run GSEA for {dsid}')
    var_cols = expset.get_varfields(dsid)
    lfc_cols = [x for x in var_cols if x.endswith('__lfc')]

    gsdict2 = get_geneset_dicts()
    for k, v in gsdict2.items():
        for kv, vv in v.items():
            print(k, kv, len(vv))
    runs = []
    lg.info(f"No lfc columns: {len(lfc_cols)}")

    for i, lfc_col in enumerate(lfc_cols):
        lg.info(f"  - processing {lfc_col}")

        rnk = expset.get_dedata_simple(dsid, lfc_col)
        rnk = rnk.set_index('gene').iloc[:, 0].sort_values()

        for j, (group_hash, gdict) in enumerate(gsdict2.items()):
            runs.append((group_hash, lfc_col, rnk, gdict))

    with Pool(threads) as P:
        results = P.map(run_one_gsea_cached, runs)

    allres = []
    
    for gh, lfc, res in sorted(results, key=lambda x: x[1]):
        res = res.rename(
            columns=dict(nes=lfc + '__nes',
                         fdr=lfc + '__fdr'))
        allres.append(
            res.melt(id_vars='set_hash', var_name='columns'))

    allres = pd.concat(allres, axis=0)
    ar2 = allres.copy()
    print(ar2.head(2).T)
    
    allres = allres.pivot(index='set_hash',
                          columns='columns',
                          values='value')

    lg.warning(f'writing to {output_file} - shape {allres.shape}')
    geneset_db = get_geneset_db()

    allres_pl = pl.from_pandas(allres.reset_index())
    allres_pl.write_parquet(output_file)


@ app.command('create_db')
def create_db():

    all_group_data = []
    all_gene_data = []

    
    geneset_folder = beehive.BASEDIR / 'geneset' / 'prep'

    lg.info(f"creating database from {geneset_folder}")
    for group_folder in geneset_folder.glob('*'):
        if not group_folder.is_dir():
            continue
        if not (group_folder / 'group.yaml').exists():
            continue

        group_info = load_yaml(group_folder / 'group.yaml')
        all_group_data.append(group_info)
        gene_data = pd.read_csv(group_folder / 'genes.tsv', sep="\t")
        all_gene_data.append(gene_data)

    all_group_data_df = pd.DataFrame(all_group_data)

    all_gene_data_df = pd.concat(all_gene_data, axis=0)

    # output...
    geneset_db = get_geneset_db()

    all_group_data_df.to_sql('groups', geneset_db, if_exists="replace")
    all_gene_data_df.to_sql('genesets', geneset_db, if_exists="replace")
    geneset_db.execute(
        """CREATE INDEX IF NOT EXISTS group_index
               ON groups (group_hash, group_title, study_title) """)

    geneset_db.execute(
        """CREATE INDEX IF NOT EXISTS genes_index
               ON genesets (geneset_hash, group_hash, title) """)

    geneset_db.commit()
    geneset_db.close()


@ app.command('import')
def import_geneset(
        base_folder: Path = typer.Argument(
            ..., file_okay=False, dir_okay=True, exists=True),
        rank_cutoff: int = typer.Option(
            250, "--rank_cutoff", "-r")):

    lg.info(f"processing {base_folder}")

    output_folder = beehive.BASEDIR / 'geneset' / 'prep'
    if not output_folder.exists():
        output_folder.mkdir(parents=True)

    for study_folder in Path(base_folder).glob('*'):
        if not study_folder.is_dir():
            continue
        if not (study_folder / 'study.yaml').exists():
            continue

        study_info = load_yaml(study_folder / 'study.yaml')
        study_hash = get_hash(**study_info)

        if 'pubmed_id' in study_info:
            study_info.update(query_pubmed(study_info['pubmed_id']))

        lg.info(f"Study: {study_folder}")

        for group_folder in study_folder.glob('*'):
            if not group_folder.is_dir():
                continue
            if not (group_folder / 'group.yaml').exists():
                continue
            lg.info(f"  - group {group_folder}")
            group_info = load_yaml(group_folder / 'group.yaml')
            group_hash = get_hash(study_hash, **group_info)
            ranktype = group_info['ranktype']
            extension = 'grp' if ranktype == 'geneset' else 'rnk'

            if 'title' not in group_info:
                group_info['title'] = " ".join(
                    group_folder.name.capitalize().split('_'))

            group_genes = {}

            group_info_2 = dict(
                study_title=study_info['title'],
                study_year=study_info['year'],
                study_author=study_info['author'],
                ranktype=group_info['ranktype'],
                organism=group_info['organism'],
                group_title=group_info['title'],
                study_hash=study_hash[:10],
                group_hash=group_hash[:10])

            group_output_folder = output_folder / group_hash[:10]
            if not group_output_folder.exists():
                group_output_folder.mkdir()

            with open(group_output_folder / 'group.yaml', 'w') as F:
                yaml.dump(group_info_2, F)

            for gsr in group_folder.glob(f'*.{extension}'):
                if ranktype == 'geneset':
                    genes = open(gsr).read().split()
                    gsr_hash = get_hash(group_hash, *genes)[:10]
                    gsr_title = " ".join(gsr.name.replace(
                        '.grp', '').capitalize().split('_'))
                    group_genes[gsr_hash] = dict(
                        genes=" ".join(genes),
                        group_hash=group_hash[:10],
                        study_hash=study_hash[:10],
                        type=ranktype,
                        direction='-',
                        title=gsr_title)
                else:
                    rnk = pd.read_csv(
                        gsr, sep="\t", header=None)\
                        .sort_values(1, ascending=False)

                    for direction in ['bottom', 'top']:

                        if direction == 'bottom':
                            genes = list(reversed(list(rnk.tail(250)[0])))
                        else:
                            genes = list(rnk.head(250)[0])

                        gsr_hash = get_hash(group_hash, *genes)[:10]
                        gsr_title = " ".join(
                            gsr.name.replace('.rnk', '')
                            .capitalize().split('_')) + ' | ' \
                            + ranktype + "-" + direction

                        group_genes[gsr_hash] = dict(
                            group_hash=group_hash[:10],
                            study_hash=study_hash[:10],
                            title=gsr_title,
                            type=ranktype,
                            direction=direction,
                            genes=" ".join(genes),
                        )

            if not group_genes:
                # no genesets found...
                continue

            group_genes_df = pd.DataFrame(group_genes).T
            group_genes_df.index.name = 'geneset_hash'
            group_genes_df = group_genes_df.reset_index()
            group_genes_df = group_genes_df[
                ['geneset_hash',  'group_hash', 'study_hash', 'title',
                 'type', 'direction', 'genes']].copy()

            lg.info(f"    - writing to {group_output_folder}", )

            group_genes_df.to_csv(group_output_folder / 'genes.tsv', sep="\t",
                                  index=False)
