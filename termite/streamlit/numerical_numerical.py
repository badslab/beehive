

import streamlit as st
import numpy as np
import pandas as pd

from termite import db
from termite.streamlit import util


def numerical_numerical(experiment, datatype, plotly_config):

    import plotly.express as px

    dimred_placeholder = st.sidebar.empty()

    candidate_genes = db.find_gene_candidates(experiment, datatype)

    exclude_num = []
    whatx, valx, datax = util.get_column(st.sidebar, 'numgene', 'X', experiment, datatype,
                                         default_gene = candidate_genes[0])
    if whatx == 'num':
        exclude_num.append(valx)

    whaty, valy, datay = util.get_column(st.sidebar, 'genenum', 'Y', experiment,
                                         datatype, exclude_num=exclude_num,
                                         default_gene = candidate_genes[1])
    if whaty == 'num':
        exclude_num.append(valy)

    whatz, valz, dataz = util.get_column(st.sidebar, 'all', 'Z', experiment,
                                         datatype, exclude_num=exclude_num,
                                         default_gene = candidate_genes[1])


    dimreds = util.find_dimred(experiment)

    if dimreds:
        def dimred_select():
            sel = st.session_state['select_dimred']
            if sel == 'Pick a DimRed':
                return
            d1, d2 = dimreds[sel]
            util.set_param('whatx', 'Numerical')
            util.set_param('whaty', 'Numerical')
            util.set_param('numx', d1)
            util.set_param('numy', d2)
            
            
        dimred_placeholder.selectbox(
            "Dim.Red", on_change=dimred_select, key='select_dimred',
            options = ["Pick a DimRed"] + list(dimreds.keys()))
        
    
    data = pd.DataFrame(
        dict(x = datax,
             y = datay,
             z = dataz))
    
    
    st.title(f'{valx} / {valy} / {valz}')

    plottype = util.selectbox_mem(
        st, "Plot type", ["scatter", "density"])

    plot_placeholder = st.empty()
    
    if plottype == 'scatter':

        with st.expander("Plot settings"):
            col1, col2 = st.columns(2)
            with col1:
                pointsize = st.slider(
                    "Point size", 1, 20, value=5)
            with col2:
                opacity = st.slider(
                    "Opacity", 0.0, 1.0, value=0.9, step=0.05)
                
        fig = px.scatter(
            data, x='x', y='y', color='z')

        fig.update_traces(
            marker=dict(size=pointsize,
                        opacity=opacity))

    elif plottype == 'density':
        fig = px.density_contour(data, x="x", y="y", color='z')

    else:
        st.write("Unknown plot type?")
        st.stop()

    fig.update_layout(
        xaxis_title=f"{whatx} {valx}",
        yaxis_title=f"{whaty} {valy}",
    )
        
    plot_placeholder.plotly_chart(fig, config=plotly_config)

