"""helpers for the gene expression app."""

from functools import partial
from genericpath import isdir, isfile
import logging
import hashlib
import copy
import sqlite3
import pandas as pd
import os
import time
from pathlib import Path
from pprint import pprint
import shutil

import typer
from typer import echo
from typing import List, Optional, Dict
import yaml

import beehive
from beehive import util, expset
from beehive.util import dict_set, query_pubmed

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


@ app.command('create_db')
def create_db(
    geneset_folder: Path = typer.Argument(
        ..., file_okay=False, dir_okay=True, exists=True), ):

    all_group_data = []
    all_gene_data = []

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
    output_db_folder = beehive.DATADIR / 'geneset_db'

    if not output_db_folder.exists():
        output_db_folder.mkdir(parents=True)

    geneset_db = sqlite3.connect(
        output_db_folder / 'geneset_db.sqlite')

    print(all_group_data_df.head(2).T)
    print(all_gene_data_df.head(2).T)

    lg.info(f"writing to: {output_db_folder / 'geneset_db.sqlite'}")
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
        output_folder: Path = typer.Argument(
            ..., file_okay=False, dir_okay=True, exists=True,
            help="output folder for prepped geneset lists"),
        rank_cutoff: int = typer.Option(
            250, "--rank_cutoff", "-r")):

    lg.info(f"processing {base_folder}")

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
