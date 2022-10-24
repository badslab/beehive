"""helpers for the gene expression app."""

from functools import partial
import hashlib
import logging
import os
import pickle
import sqlite3

import pandas as pd
from pathlib import Path
from pprint import pprint
import gseapy as gp
import polars as pl
import yaml


import typer
from typer import echo
from typing import List, Optional, Dict

import beehive
from beehive import util, expset
from beehive.util import dict_set, diskcache, query_pubmed, \
    get_geneset_db
from beehive import util, expset


app = typer.Typer()

lg = logging.getLogger(__name__)
lg.setLevel(logging.INFO)


@app.command('ds')
def ds():
    datasets = expset.get_datasets()
    for name, ds in datasets.items():
        print("\t".join([
            name, ds['short_author'], ds['short_title'],
            ds['datatype'], ]))


@app.command('de')
def de(dsid: str,
       ):
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

    gd = expset.get_gsea_data(dsid, decol)

    del gd['geneset_hash']
    del gd['study_hash']
    del gd['set_hash']

    print(gd['de study_author title fdr nes'.split()].head(20))
