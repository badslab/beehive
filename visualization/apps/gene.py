
from functools import partial

import streamlit as st
import pandas as pd

from termite.vis import util as vutil
import blocks, views
from termite.vis import data as vdata
    
def gene():

    
    dataset_rec = blocks.get_dataset()
    gene = blocks.get_gene("gene", dataset_rec['dataset_id'])
    expr = vdata.get_expr(dataset_rec['dataset_id'], gene)

    views.Multi(
        views.Ecdf(data=expr, subject=gene),
        views.Histogram(data=expr, subject=gene)
        ).view()

    
def gene2():

    dataset_rec = blocks.get_dataset()
    ds_id = dataset_rec['dataset_id']
    gene1 = blocks.get_gene("gene", dataset_id=ds_id)
    gene2 = blocks.get_gene("gene2", dataset_id=ds_id)
    expr1 = vdata.get_expr(ds_id, gene1)
    expr2 = vdata.get_expr(ds_id, gene2)

    data = pd.DataFrame({gene1:expr1, gene2:expr2})

    views.Multi(
        views.DensityHeatmap(data=data, subject=f"{gene1} vs {gene2}"),
        views.DensityHeatmap2(data=data, subject=f"{gene1} vs {gene2}")
    ).view()

