

from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator


from termite import db


@st.cache_data
def df_to_tsv(df):
    rv = df.to_csv(sep="\t").encode('utf-8')
    return rv


def set_param(key, val):
    qp = st.experimental_get_query_params()
    qp[key] = val
    st.session_state[key] = val

    
def suggest_genes(
        gene: str,
        exp_id: int,
        inject_in: str = 'gene') -> None:
    
    # probably gene not found: suggest a few:
    candidate_genes = list(db.fuzzy_gene_search(exp_id, gene))
    candidate_genes = ["Please pick one"] + candidate_genes
    
    def _gene_select():
        if st.session_state['suggest_a_gene'] != "Please pick one":
            picked = st.session_state['suggest_a_gene']
            qp = st.experimental_get_query_params()
            qp[inject_in] = picked
            st.experimental_set_query_params(**qp)
            st.session_state[inject_in] = picked
        
    st.selectbox(
        f'Gene "{gene}" not found, did you mean:',
        candidate_genes, key='suggest_a_gene',
        on_change=_gene_select)

    
def selectbox_mem(
        context: DeltaGenerator,
        label: str,
        options: list[str],
        default: Optional[str] = None,
        index: int =0,
        key: Optional[str] = None) -> str:
    
    if key is None:
        key = label.lower().replace(' ', '_')
        
    options = list(options)
    qp = st.experimental_get_query_params()

    def update_query_param():
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    if default is not None and default in options:
        idx = options.index(default)
    else:
        idx = index
        
    if key in qp:
        qval = qp[key][0]
        if qval in options:
            idx = options.index(qval)
        
    return context.selectbox(label, options, key=key, index=idx,
                             on_change=update_query_param)


def textbox_mem(
        context: DeltaGenerator,
        label: str,
        default: str="",
        key: Optional[str]=None) -> str:

    if key is None:
        key = label.lower().replace(' ', '_')
        
    qp = st.experimental_get_query_params()

    def update_query_param():
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    # set default value via session state - otherwise we get
    # errors when we update this through other methods
    if key in qp:
        dvalue = qp[key][0]
    else:
        dvalue = default
        
    return context.text_input(label, key=key, value=dvalue,
                             on_change=update_query_param)


def get_column(
        context: DeltaGenerator,
        which_types: str,
        label_suffix: str,
        exp_id: int,
        default_gene: Optional[str] = None,
        default_num: Optional[str] = None,
        default_cat: Optional[str] = None,
        key_suffix: Optional[str] = None,
        exclude_num: Optional[List[str]] = None,
        exclude_cat: Optional[List[str]] = None) \
        -> Tuple[str, str, pd.Series]:
    
    """ Generic function to get a column of data from the
        database. """


    # prepare numnames
    all_numnames = db.get_obs_num_names(exp_id)
    if exclude_num:
        all_numnames = set(all_numnames) - set(exclude_num)
    all_numnames = list(sorted(all_numnames))

    # prepare catnames
    all_catnames = db.get_obs_cat_names(exp_id)
    if exclude_cat:
        all_catnames = set(all_catnames) - set(exclude_cat)
    all_catnames = list(sorted(all_catnames))

    if key_suffix is None:
        key_suffix = "".join(label_suffix.split()).lower()


    if which_types in ['genenum', 'num']:
        what_types = ['Gene', 'Gene/log1p',
                      'Numerical', 'Numerical/log1p']
    elif which_types == 'numgene':
        what_types = ['Numerical', 'Numerical/log1p',
                      'Gene', 'Gene/log1p']              
    elif which_types == 'gene':
        what_types = ['Gene', 'Gene/log1p']
    elif which_types == 'cat':
        what_types = ['Categorical']
    else:
        what_types = ['Gene', 'Gene/log1p',
                      'Numerical', 'Numerical/log1p',
                      'Categorical']
    
    if len(what_types) > 1:
        col1, col2 = context.columns(2)
        what = selectbox_mem(
            context=col1,
            label="Plot " + label_suffix,
            options=what_types,
            key='what' + key_suffix)
    else:
        what = what_types[0]
        col2 = context

        
    def transform(_what, _data):
        if 'log1p' in _what:
            _data = np.log1p(_data)
        elif 'sqrt' in _what:
            _data = np.sqrt(_data)
        return _data
    
    if "Gene" in what:
        gene = textbox_mem(
            context=col2,
            label="Gene " + label_suffix,
            key="gene" + key_suffix,
            default=default_gene)
        
        data = db.get_expr_gene(exp_id, gene)
        if len(data) == 0:
            # gene not found - suggest a few candidateas
            suggest_genes(gene, exp_id, inject_in='gene' + key_suffix)
            st.stop()
        data = transform(what, data)
        return "gene", gene, data

    elif 'Numerical' in what:
        num = selectbox_mem(
            context=col2,
            label="Numerical " + label_suffix,
            options=all_numnames,
            key="num" + key_suffix)
        data = db.get_obs_num(exp_id, num)
        data = transform(what, data)
        return "num", num, data

    else:
        cat = selectbox_mem(
            context=col2,
            label="Categorical " + label_suffix,
            options=all_catnames,
            key="cat" + key_suffix)
        data = db.get_obs_cat(exp_id, cat)
        return "cat", cat, data

