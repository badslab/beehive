"""helpers for the gene expression app."""

import copy
from datetime import datetime

import hashlib
import logging
import pickle
from pathlib import Path

import pandas as pd
import polars as pl
import typer
import yaml

import beehive
from beehive import util
from beehive.util import get_geneset_db, query_pubmed

app = typer.Typer()

lg = logging.getLogger(__name__)
lg.setLevel(logging.DEBUG)


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


def get_geneset_groups(dsid, organism):
    gsdb = get_geneset_db(dsid)
    groups = pd.read_sql(
        f'''SELECT *
             FROM groups
            WHERE organism="{organism}"
        ''', gsdb)
    yield from groups.iterrows()


def get_geneset_genes(dsid, group_hash):
    gsdb = get_geneset_db(dsid)
    gsets = pd.read_sql(
        f'''SELECT * FROM genesets
            WHERE group_hash = '{group_hash}'
        ''', gsdb)
    return gsets


def run_one_gsea(args):

    import logging

    import gseapy as gp
    logging.getLogger("gseapy").setLevel(logging.ERROR)

    group_hash, lfc_col, rnk, genedict = args

    results = gp.prerank(
        rnk=rnk, gene_sets=genedict,
        threads=3, min_size=5,
        max_size=1000, permutation_num=5000,
        outdir=None, seed=42,
        verbose=True)

    print(results.res2d)
    res2 = results.res2d[['Term', 'NES', 'FDR q-val', 'Lead_genes']]
    res2.columns = ['geneset_hash', 'nes', 'fdr', 'lead_genes']
    res2['lead_genes'] = res2['lead_genes'].str.replace(';', ' ')
    res2.set_index('geneset_hash')
    lg.warning(f"Finished one gsea {group_hash} - {lfc_col}")
    return group_hash, lfc_col, res2


def run_one_gsea_cached(args):
    uid = util.UID(*args, length=12)
    cache_folder = util.get_geneset_folder() / 'cache' / 'gsea'

    if not cache_folder.exists():
        try:
            cache_folder.mkdir(parents=True)
        except FileExistsError:
            # in a multi threaded run -another thread might
            # attempt the same - let's not crash..
            pass

    cache_file = cache_folder / uid

    if cache_file.exists():
        with open(cache_file, 'rb') as F:
            lg.debug("return from cache")
            rv = pickle.load(F)
    else:
        rv = run_one_gsea(args)
        lg.debug(f"saving to cache: {uid}")
        with open(cache_file, 'wb') as F:
            pickle.dump(rv, F)

    return rv


def get_geneset_dicts(dsid: str,
                      organism: str) -> dict:
    """Return a dict of dictionaries of a given
       organism
    """

    rv = {}

    for _, g in get_geneset_groups(dsid, organism=organism):

        gsets = get_geneset_genes(dsid, g['group_hash'])
        lg.info(
            f"    | {g['study_title'][:40]} "
            f" | {g['group_title'][:40]}")

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

    from beehive import expset

    dsyaml = expset.get_dataset(dsid, None)
    organism = dsyaml['organism']
    lg.info(f"Organism {organism}")

    lg.info(f'Run GSEA for {dsid}')
    var_cols = expset.get_varfields(dsid)
    lfc_cols = [x for x in var_cols if x.endswith('__lfc')]

    gsdict2 = get_geneset_dicts(dsid, organism=organism)
    runs = []
    lg.info(f"No lfc columns: {len(lfc_cols)}")

    for i, lfc_col in enumerate(lfc_cols):

        rnk = expset.get_dedata_simple(dsid, lfc_col)
        lg.info(f"  - processing {lfc_col} ({len(rnk)} genes)")
        rnk = rnk.set_index('gene').iloc[:, 0].sort_values()

        for j, (group_hash, gdict) in enumerate(gsdict2.items()):
            runs.append((group_hash, lfc_col, rnk, gdict))

    with Pool(threads) as P:
        results = P.map(run_one_gsea_cached, runs)

    allres_raw = []
    for gh, lfc, res in sorted(results, key=lambda x: x[1]):
        res['column'] = lfc
        allres_raw.append(res)

    allres = pd.concat(allres_raw, axis=0)
    geneset_db = get_geneset_db(dsid)
    allres.to_sql('gsea', geneset_db,
                  if_exists="replace", index=False)

    geneset_db.execute(
        """CREATE INDEX IF NOT EXISTS gsea_index
               ON gsea (geneset_hash, column, nes, fdr) """)

    geneset_db.commit()
    geneset_db.close()


@ app.command('create-db')
def create_db(
        dsid: str = typer.Argument(..., help='Dataset'),):

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
    geneset_db = get_geneset_db(dsid)

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


@app.command('gsea-export')
def gsea_export(
        dsid: str = typer.Argument(..., help='Dataset'),):

    output_path = util.get_geneset_folder() / "gsea" / dsid
    if not output_path.exists():
        output_path.mkdir(parents=True)

    gsdb = get_geneset_db(dsid)

    columns = pd.read_sql(
        '''select distinct(column) from gsea''', gsdb)
    for _, column in columns['column'].iteritems():
        lg.info(f"exporting {column}")
        sql = \
            f''' SELECT
                    gsea.geneset_hash,
                    gsea.nes,
                    gsea.fdr,
                    gsea.lead_genes,
                    gsea.column as de_column,
                    gs.title as geneset_title,
                    gs.type as geneset_type,
                    gs.direction as direction,
                    gs.genes as genes,
                    gr.organism as organism,
                    gr.study_title as study_title,
                    gr.study_author as study_author,
                    gr.study_year as study_year
                FROM genesets as gs, groups as gr,
                     gsea as gsea
                WHERE gs.group_hash = gr.group_hash
                  AND gsea.geneset_hash = gs.geneset_hash
                  AND gsea.column = "{column}"
                ORDER BY gsea.fdr
            '''
        data = pd.read_sql(sql, gsdb)
        data['no_genes'] = data['genes'].str.split().apply(len)
        data = data.set_index('geneset_hash')
        assert data.index.is_unique

        output_name = f"gsea__{dsid}__{column}"
        output_csv_file = output_path / f"{output_name}.tsv"
        output_xlsx_file = output_path / f"{output_name}.xlsx"
        print(f'export to {output_name}')
        data.to_csv(output_csv_file, sep="\t")
        data.to_excel(output_xlsx_file)



@diskcache()
def import_one_from_enrichr(organism, geneset):
    lg.info("downloading from enrichr")
    return gp.get_library(name=geneset, organism=organism)
    
@ app.command('import-enrichr')
def import_enrichr(
        geneset: str,
        organism = typer.Option('Human', '--organism', '-o'),
        ):
    
    """ @Todo: organism support? Now sticking with human
    """
    print(organism, geneset)
    enrgenes = import_one_from_enrichr(organism, geneset)

    output_folder = beehive.BASEDIR / 'geneset' / 'prep'
    if not output_folder.exists():
        output_folder.mkdir(parents=True)

    now = datetime.now()
    study_info = dict(title='enrichr',
                      year=now.year,
                      author='enrichr')
    study_hash = get_hash(**study_info)
    group_info = dict(title=geneset,
                      organism=organism.lower(),
                      ranktype='geneset',
                      year=now.year,
                      author='enrichr')
    group_hash = get_hash(study_hash, **group_info)
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

    group_genes = {}
    
    for gsr_title, genes in enrgenes.items():
        gsr_hash = get_hash(group_hash, *genes)[:10]
        group_genes[gsr_hash] = dict(
            genes=" ".join(genes),
            group_hash=group_hash[:10],
            study_hash=study_hash[:10],
            type='geneset', direction='-',
            title=gsr_title)

    group_genes_df = pd.DataFrame(group_genes).T
    group_genes_df.index.name = 'geneset_hash'
    group_genes_df = group_genes_df.reset_index()
    group_genes_df = group_genes_df[
        ['geneset_hash', 'group_hash', 'study_hash', 'title',
         'type', 'direction', 'genes']].copy()

    lg.info(f"    - writing to {group_output_folder}", )

    group_genes_df.to_csv(group_output_folder / 'genes.tsv', sep="\t",
                          index=False)


    
@ app.command('import')
def import_geneset(
        rank_cutoff: int = typer.Option(
            250, "--rank_cutoff", "-r")):

    base_folder = util.get_geneset_folder() / "db"
    if not base_folder.exists():
        lg.error("Can not find geneset db: {base_folder}")
        exit()
    lg.info(f"Processing {base_folder}")

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
                ['geneset_hash', 'group_hash', 'study_hash', 'title',
                 'type', 'direction', 'genes']].copy()

            lg.info(f"    - writing to {group_output_folder}", )

            group_genes_df.to_csv(group_output_folder / 'genes.tsv', sep="\t",
                                  index=False)
