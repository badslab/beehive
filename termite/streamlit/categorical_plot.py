

import streamlit as st
import numpy as np


from termite import db
from termite.streamlit import util


def categorical_plot(experiment, datatype, plotly_config):
    
    import plotly.express as px

    
    all_catnames = db.get_obs_cat_names(experiment)
    catname = util.selectbox_mem(
        st.sidebar, "Categorical metadata", all_catnames,
        key='cat1')

    # what to plot as numerical - a gene, or numerical metadata
    gene_or_num = util.selectbox_mem(
        st.sidebar, 'Gene or Numerical', ['Gene', 'Numerical metadata'])

    if gene_or_num == 'Gene':
        gene = util.textbox_mem(st.sidebar, "Gene", default='APOE', key='gene1')
        st.title(f'Plot Gene "{gene}" vs "{catname}"')
        data = db.get_expr_obscat(experiment, datatype, gene, catname)
        num_var_name = f'{gene} expression ({datatype})'
        file_name = f"cateplot_{experiment}_{datatype}_{gene}_{catname}.tsv"
    else:
        all_numnames = list(sorted(db.get_obs_num_names(experiment)))
        numname = st.sidebar.selectbox(
            "Numerical metadata", all_numnames, key='num1')
        st.title(f'Plot Variable "{numname}" vs "{catname}"')
        data = db.get_numvar_obscat(experiment, numname, catname)
        num_var_name = numname
        file_name = f"catplot_{experiment}_{datatype}_{numname}_{catname}.tsv"

        
    data = data.sort_values(by='cat')
    
    if gene_or_num == 'Gene' and len(data) == 0:
        util.suggest_genes(gene, experiment, datatype, 'gene1')
        st.stop()

    plottype = st.sidebar.selectbox(
        'Plot type', ['Violin', 'Box', 'Reversed ecdf', 'Bar'])

    if plottype == 'Violin':
        fig = px.violin(data, x='cat', y='num', box=True, points='outliers')
        fig.update_layout(                
            xaxis_title=f"{catname}",
            yaxis_title=f"{num_var_name}",
        )

    elif plottype == 'Box':
        fig = px.box(data, x='cat', y='num')
        fig.update_layout(                
            xaxis_title=f"{catname}",
            yaxis_title=f"{num_var_name}")
    elif plottype == 'Reversed ecdf':
        fig = px.ecdf(data, y='num', orientation='h',
                      color='cat')
        fig.update_layout(                
            yaxis_title=f"{num_var_name}")

    elif plottype == 'Bar':
        from scipy.stats import sem
        col1, col2 = st.columns(2)
        with col1:
            stat = st.selectbox(
                'Aggregate statistic', ['Mean', 'Median'])
        with col2:
            errorbar = st.selectbox(
                'Error bar', ['StdError', 'Stdev', 'None'])
        
        agg = data.groupby('cat')['num']\
                    .agg(Mean=np.mean,
                         Median=np.median,
                         Stdev=np.std,
                         StdError=sem,)\
                    .reset_index()
        
        if errorbar != 'None':
            fig = px.bar(agg, x='cat', y=stat,
                         error_y=errorbar)
        else:
            fig = px.bar(agg, x='cat', y=stat)

        fig.update_layout(                
            xaxis_title=f"{catname}",
            yaxis_title=f"{stat} {num_var_name}",
        )


    else:
        st.write("Unknow plot type?")
        st.stop()

        
    st.plotly_chart(fig, config=plotly_config)
        
    st.download_button(
        "Download data",
        util.df_to_tsv(data),
        file_name=file_name,
        mime="text/tsv")
                                    
