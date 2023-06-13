

import streamlit as st
import numpy as np

from termite import db
from termite.streamlit import util


def gene_vs_category(experiment, datatype, plotly_config):
    
    import plotly.express as px
    
    gene = st.sidebar.text_input("Gene", value='APOE', key='gene')


    all_catnames = db.get_obs_cat_names(experiment)
    catname = st.sidebar.selectbox(
        "Categorical metadata", all_catnames)
    
    st.title(f"Plot {gene} vs {catname}")
    
    counts = db.get_expr_obscat(experiment, datatype, gene, catname)

    counts = counts.sort_values(by='cat')
    
    if len(counts) == 0:
        # probably gene not found: suggest a few:
        candidate_genes = db.fuzzy_gene_search(experiment, gene)
        def _gene_select():
            st.session_state['gene'] = st.session_state['new_gene']
        st.selectbox(
            'Gene not found, did you mean:',
            candidate_genes, key='new_gene',
            on_change=_gene_select)

        st.stop()

    plottype = st.sidebar.selectbox(
        'Plot type', ['Violin', 'Box', 'Reversed ecdf', 'Bar'])

    if plottype == 'Violin':
        fig = px.violin(counts, x='cat', y='expr', box=True, points='outliers')
    elif plottype == 'Box':
        fig = px.box(counts, x='cat', y='expr')
    elif plottype == 'Reversed ecdf':
        fig = px.ecdf(counts, y='expr', orientation='h',
                      color='cat')
    elif plottype == 'Bar':
        from scipy.stats import sem
        col1, col2 = st.columns(2)
        with col1:
            stat = st.selectbox(
                'Aggregate statistic', ['Mean', 'Median'])
        with col2:
            errorbar = st.selectbox(
                'Error bar', ['Stdev', 'StdError'])
        
        agg = counts.groupby('cat')['expr']\
                    .agg(Mean=np.mean,
                         Median=np.median,
                         Stdev=np.std,
                         StdError=sem,)\
                    .reset_index()

        fig = px.bar(agg, x='cat', y=stat,
                     error_y=errorbar)

    else:
        st.stop()
        
    st.plotly_chart(fig, config=plotly_config)
        
    st.download_button(
        "Download data",
        util.df_to_tsv(counts),
        file_name=f"categ_gene_{experiment}_{datatype}_{gene}_{catname}.tsv",
        mime="text/tsv")
                                    
