
from functools import partial

import streamlit as st
import pandas as pd

from termite.vis import util as vutil
import blocks, views
from termite.vis import data as vdata
    
def gene():

    
    dataset_rec = blocks.get_dataset()
    gene, expr = blocks.get_expr("gene", dataset_rec['dataset_id'])

    views.Multi(
        views.SnsEcdf(data=expr, subject=gene),
        views.Histogram(data=expr, subject=gene)
        ).view()

    
def gene2():

    dataset_rec = blocks.get_dataset()
    exp_id=dataset_rec['experiment_id']
    ds_id=dataset_rec['dataset_id']
    gene1, expr1 = blocks.get_expr("gene", dataset_id=ds_id)
    gene2, expr2 = blocks.get_expr("gene2", dataset_id=ds_id)

    data = pd.DataFrame({gene1:expr1, gene2:expr2})

    views.Multi(
        views.DensityHeatmap(data=data, subject=f"{gene1} vs {gene2}"),
        views.DensityHeatmap2(data=data, subject=f"{gene1} vs {gene2}")
    ).view()

