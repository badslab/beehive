

import streamlit as st
import numpy as np


from termite import db
from termite.streamlit import util


def categorical_plot(experiment, datatype, plotly_config):

    
    import pandas as pd
    import plotly.express as px
    from scipy.stats import sem


    _, catname, catcol = util.get_column(
        st.sidebar, 'cat', 'X',
        experiment=experiment, datatype=datatype)
    
    typey, coly, datay = util.get_column(
        st.sidebar, 'num', 'Y', 
        experiment=experiment, datatype=datatype)

    data = pd.DataFrame(dict(
        cat = catcol,
        num = datay))

    agg = data.groupby('cat')['num']\
              .agg(Mean=np.mean,
                   Median=np.median,
                   Stdev=np.std,
                   StdError=sem,)\
              .reset_index()

    filename_raw = f"catplot_{experiment}_{catname}_{coly}_raw.tsv"
    filename_agg = f"catplot_{experiment}_{catname}_{coly}_agg.tsv"
    
    plottype = st.sidebar.selectbox(
        'Plot type', ['Violin', 'Box', 'Reversed ecdf', 'Bar'])

    if plottype == 'Violin':
        fig = px.violin(data, x='cat', y='num', box=True, points='outliers')
        fig.update_layout(                
            xaxis_title=f"{catname}",
            yaxis_title=f"{coly}")

    elif plottype == 'Box':
        fig = px.box(data, x='cat', y='num')
        fig.update_layout(                
            xaxis_title=f"{catname}",
            yaxis_title=f"{coly}")
    elif plottype == 'Reversed ecdf':
        fig = px.ecdf(data, y='num', orientation='h',
                      color='cat')
        fig.update_layout(                
            yaxis_title=f"{coly}")

    elif plottype == 'Bar':
        col1, col2 = st.columns(2)
        with col1:
            stat = st.selectbox(
                'Aggregate statistic', ['Mean', 'Median'])
        with col2:
            errorbar = st.selectbox(
                'Error bar', ['StdError', 'Stdev', 'None'])
        
        
        if errorbar != 'None':
            fig = px.bar(agg, x='cat', y=stat,
                         error_y=errorbar)
        else:
            fig = px.bar(agg, x='cat', y=stat)

        fig.update_layout(                
            xaxis_title=f"{catname}",
            yaxis_title=f"{stat} {coly}",
        )


    else:
        st.write("Unknow plot type?")
        st.stop()

        
    st.plotly_chart(fig, config=plotly_config)

    
    with st.empty():
        if st.button('Prepare data downloads'):
            with st.container():
                st.download_button(
                    "Download raw data",
                    util.df_to_tsv(data),
                    file_name=filename_raw,
                    mime="text/tsv")

                st.download_button(
                    "Download aggregated data",
                    util.df_to_tsv(agg),
                    file_name=filename_agg,
                    mime="text/tsv")
                                    
