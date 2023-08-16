
from functools import partial

import pandas as pd
import streamlit as st
import views

from termite.vis import util as vutil
from termite import util

def dimred():

    import blocks
    import views

    dataset_rec = blocks.get_dataset()
    ds_id =  dataset_rec['dataset_id']
    exp_id = dataset_rec['experiment_id']
    dimred = blocks.get_dimred_name(exp_id)

    color_what = vutil.selectbox_mem(
        st.sidebar, "Color on",
        options=["Gene", "Numeric md", "Categorical md"],
        key='colon')

    #if color_what == 'gene'

    gene = blocks.get_gene('gene', dataset_id =ds_id)
    gene_expr = blocks.get_expr(dataset_id=ds_id, gene=gene)
    dim1 = blocks.get_obs_num(exp_id=exp_id, colname=f"{dimred}/00")
    dim2 = blocks.get_obs_num(exp_id=exp_id, colname=f"{dimred}/01")

    data = pd.DataFrame(dict(dim1=dim1, dim2=dim2, zcol=gene_expr))

    
    scatter = partial(
        views.scatter,
        labels=dict(
            dim1=f"{dimred}/00", 
            dim2=f"{dimred}/01",
            zcol=gene)
        )

    view = vutil.DictSelect(
        "View",
        dict(Scatter=scatter,
             ))


    ctxPlot = st.empty()
    ctxParam = st.container()
    view(data, ctxPlot, ctxParam)
    
