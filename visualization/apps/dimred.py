
from functools import partial

import pandas as pd
import streamlit as st
import views

from termite.vis import util

def dimred():

    import blocks
    import views

    dataset_rec = blocks.get_dataset()
    ds_id =  dataset_rec['dataset_id']
    exp_id = dataset_rec['experiment_id']
    dimred = blocks.get_dimred_name(exp_id)


    col_what, col_name, col = blocks.get_gcn(
        "Color on", ds_id, exp_id, key='color')

    filter_what, filter_name, filter_col = blocks.get_gcn(
        "Filter on", ds_id, exp_id, key='filter',
        allow_nothing=True)

    
    dim1 = blocks.get_obs_num(exp_id=exp_id, colname=f"{dimred}/00")
    dim2 = blocks.get_obs_num(exp_id=exp_id, colname=f"{dimred}/01")

    
    data = pd.DataFrame(
        {'dim1': dim1,
         'dim2': dim2,
         })
    
    data = data.merge(col, how='inner', left_index=True, right_index=True)

    if filter_what != "None":
        data = data.merge(filter_col, how='inner', left_index=True, right_index=True)
        
    views.Multi(
        views.Scatter(data=data, xloggable=False, yloggable=False,
                      zloggable=(col_what != 'Cat'),
                      pointsize=3, subject=f"{dimred}"),

    ).view()
    
