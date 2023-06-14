

import streamlit as st
import numpy as np
import pandas as pd

from termite import db
from termite.streamlit import util


def numerical_numerical(experiment, datatype, plotly_config):

    import plotly.express as px

    exclude_num = []
    whatx, valx, datax = util.get_column(st.sidebar, 'num', 'X', experiment, datatype)
    if whatx == 'num':
        exclude_num.append(valx)
        
    whaty, valy, datay = util.get_column(st.sidebar, 'num', 'Y', experiment,
                                         datatype, exclude_num=exclude_num)
    if whaty == 'num':
        exclude_num.append(valy)

    whatz, valz, dataz = util.get_column(st.sidebar, 'all', 'Z', experiment,
                                         datatype, exclude_num=exclude_num)

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

