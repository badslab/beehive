

import streamlit as st
from termite import db


@st.cache_data
def df_to_tsv(df):
    return df.to_csv(sep="\t").encode('utf-8')


def suggest_genes(gene, experiment, datatype):
    # probably gene not found: suggest a few:
    candidate_genes = db.fuzzy_gene_search(experiment, datatype, gene)
    def _gene_select():
        st.session_state['gene'] = st.session_state['new_gene']
        
    st.selectbox(
        'Gene not found, did you mean:',
        candidate_genes, key='new_gene',
        on_change=_gene_select)


def selectbox_mem(context, label, options, index=0):

    key = label.lower().replace(' ', '_')
    options = list(options)
    qp = st.experimental_get_query_params()

    def update_query_param():
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    idx = index
    if key in qp:
        qval = qp[key][0]
        if qval in options:
            idx = options.index(qval)
        
    return context.selectbox(label, options, key=key, index=idx,
                             on_change=update_query_param)
    

    
