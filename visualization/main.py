"""Streamlit visualization for termite."""

import logging

import apps.dataset
import apps.dimred
import apps.gene
import apps.gentab
import apps.helper
import streamlit as st
from rich.logging import RichHandler

from termite.vis import util

FORMAT = "%(message)s"
logging.basicConfig(
    level=0,
    format=FORMAT,
    datefmt="[%X]",
    handlers=[RichHandler()])

lg = logging.getLogger(__name__)
lg.info("Start Termite Visualization Streamlit App!")

st.set_page_config(layout="wide")

# Make  captions look like labels - for layout purposes.
st.sidebar.markdown("""<style>

div[data-testid="stCheckbox"] > label{
  text-align: end;
  align-items: center;
  justify-content: right;
}

[data-testid="stHorizontalBlock"] {
  gap: 0.2rem;
}

[data-testid="stCaptionContainer"] p {
  font-family: "Source Sans Pro", sans-serif;
  font-size: 14px;
  padding: 0px;
  margin-bottom: 0.25rem;
  font-weight: 400;
  line-height: 24px;
  color: rgb(49, 51, 63); }
</style>



""", unsafe_allow_html=True)

all_apps = {
    "dataset": ("Dataset",             apps.dataset.dataset),
    "general": ("Generic Table Maker", apps.gentab.gentab),
    "categ":   ("Categorical MD",      apps.dataset.categorical),
    "cat2":    ("Two Categorical MD",  apps.dataset.catcat),
    "num2":    ("Two Numerical MD",    apps.dataset.numnum),
    "numer":   ("Numerical MD",        apps.dataset.numerical),
    "gene":    ("Gene/Exp",            apps.gene.gene),
    "gene2":   ("Two Genes",           apps.gene.gene2),
    "dimred":  ("Dimred",              apps.dimred.dimred),
    "tblcnt":  ("Table count",         apps.helper.tablepeek),
    "sql":     ("Raw SQL",             apps.helper.raw_sql),
}


col1, col2 = st.sidebar.columns([3, 7])
col1.caption('**Application**')
appname, app = util.DictSelect(
    'Application', all_apps, key='app',
    label_visibility='collapsed', context=col2,
    format_func=lambda x: all_apps[x][0])

app()
