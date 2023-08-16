
import logging

from rich.logging import RichHandler
import streamlit as st

from termite.vis import util

import apps.dimred
import apps.dataset
import apps.gene
import apps.helper


FORMAT = "%(message)s"
logging.basicConfig(
    level=0,
    format=FORMAT,
    datefmt="[%X]",
    handlers=[RichHandler()])

lg = logging.getLogger(__name__)
lg.info("Start Termite Visualization Streamlit App!")

st.set_page_config(layout="wide")

all_apps = {
    "dataset": ("Dataset",        apps.dataset.dataset),
    "categ":   ("Categorical MD", apps.dataset.categorical),
    "cat2":    ("Two Categorical MD", apps.dataset.catcat),
    "numer":   ("Numerical MD",   apps.dataset.numerical),
    "gene":    ("Gene/Exp",       apps.gene.gene),
    "gene2":   ("Two Genes",      apps.gene.gene2),
    "dimred":  ("Dimred",         apps.dimred.dimred),
    "tblcnt":  ("Table count",    apps.helper.tablepeek),
    "sql":     ("Raw SQL",        apps.helper.raw_sql),
}


appname, app = util.DictSelect(
    'Application', all_apps, key='app',
    format_func=lambda x: all_apps[x][0])

app()
    
