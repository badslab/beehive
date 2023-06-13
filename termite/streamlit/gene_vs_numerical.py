

import streamlit as st
import numpy as np

from termite import db
from termite.streamlit import util


def gene_vs_numerical(experiment, datatype, plotly_config):

    import plotly.express as px

    all_numnames = list(sorted(db.get_obs_num_names(experiment)))

    gene = st.sidebar.text_input("Gene", value='APOE', key='gene')

    expr_scale = util.selectbox_mem(
        st.sidebar, "Expression scale", ['Unscaled', 'log1p'])
        
    numname = st.sidebar.selectbox(
        "Numerical metadata", all_numnames, key='cat1')
    
    num_scale = util.selectbox_mem(
        st.sidebar, "Num. MD scale", ['Unscaled', 'log1p'])

    st.title(f'"{gene}" vs "{numname}"')

    counts = db.get_expr_obsnum(
        experiment, datatype, gene, numname)
    counts = counts.sort_values(by=['num1'])

    ylabel_extra = ""
    if expr_scale == 'log1p':
        counts['expr'] = np.log1p(counts['expr'])
        ylabel_extra = "(log1p)"
        
    xlabel_extra = ""
    if num_scale == 'log1p':
        counts['num1'] = np.log1p(counts['num1'])
        xlabel_extra = "(log1p)"
        
    if len(counts) == 0:
        util.suggest_genes(gene, experiment, datatype)
        st.stop()

    plottype = util.selectbox_mem(
        st, "Plot type", ["scatter", "density"])


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
            counts, x='num1', y='expr')

        fig.update_traces(
            marker=dict(size=pointsize,
                        opacity=opacity))

    elif plottype == 'density':
        fig = px.density_contour(counts, x="num1", y="expr")
        
    else:
        st.write("Unknown plot type?")
        st.stop()

    fig.update_layout(
        xaxis_title=f"{numname} {xlabel_extra}",
        yaxis_title=f"{gene} expression {ylabel_extra}",
    )
        
    st.plotly_chart(fig, config=plotly_config)

