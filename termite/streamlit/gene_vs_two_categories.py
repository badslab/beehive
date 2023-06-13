

import streamlit as st
import numpy as np

from termite import db
from termite.streamlit import util


def gene_vs_two_categories(experiment, datatype, plotly_config):

    import plotly.express as px
    from scipy.stats import sem

    gene = st.sidebar.text_input("Gene", value='APOE', key='gene')

    all_catnames = list(sorted(db.get_obs_cat_names(experiment)))
    catname1 = st.sidebar.selectbox(
        "Categorical metadata 1", all_catnames, key='cat1')

    all_catnames_2 = list(sorted(set(all_catnames) - set([catname1])))
    catname2 = st.sidebar.selectbox(
        "Categorical metadata 2", all_catnames_2, key='cat2')
    

    
    st.title(f'"{gene}" vs "{catname1}" and "{catname2}"')

    
    counts = db.get_expr_obscat_two(experiment, datatype, gene, catname1, catname2)
    counts = counts.sort_values(by=['cat1', 'cat2'])

    
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

        
    agg = counts.groupby(['cat1', 'cat2'])['expr']\
            .agg(Mean=np.mean,
                 Median=np.median,
                 Stdev=np.std,
                 StdError=sem,)\
            .reset_index()

    
    fig = px.bar(agg, x='cat1', y='Mean', color='cat2',
                 barmode='group', error_y='StdError')

        
    st.plotly_chart(fig, config=plotly_config)

