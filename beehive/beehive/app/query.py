"""helpers for the gene expression app."""

import logging

import typer

app = typer.Typer()

lg = logging.getLogger(__name__)
lg.setLevel(logging.INFO)


@app.command('ds')
def ds():
    from beehive import expset
    datasets = expset.get_datasets()
    for name, ds in datasets.items():
        print("\t".join([
            name, ds['short_author'], ds['short_title'],
            ds['datatype'], ]))


@app.command('de')
def de(dsid: str,
       ):
    from beehive import expset
    cols = expset.get_varfields(dsid)
    if 'gene' in cols:
        cols.remove('gene')
    else:
        cols = cols[:-1]

    for c in cols:
        print(c)


@app.command('gsea')
def gsea(dsid: str,
         decol: str = typer.Argument(None),):
    from beehive import expset
    gd = expset.get_gsea_data(dsid, decol)

    del gd['geneset_hash']
    del gd['study_hash']
    del gd['set_hash']

    print(gd['de study_author title fdr nes'.split()].head(20))
