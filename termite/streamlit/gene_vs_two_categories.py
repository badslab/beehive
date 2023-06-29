

import streamlit as st
import numpy as np


from termite.streamlit import util
from termite import db

def gene_vs_two_categories(
        experiment: str,
        exp_id: int,
        plotly_config: dict):

    import plotly.express as px
    import pandas as pd
    from scipy.stats import sem

    candidate_genes = db.find_gene_candidates(exp_id)

    xtype, num, numdata = util.get_column(
        st.sidebar, "num", "x", exp_id,
        default_gene=candidate_genes[0],)
    
    _, cat1name, cat1data = util.get_column(
        st.sidebar, "cat", "1", exp_id)
    _, cat2name, cat2data = util.get_column(
        st.sidebar, "cat", "2", exp_id,
        exclude_cat = [cat1name] )
    
    
    filename_raw = f"cat2plot_{experiment}_{cat1name}_{cat2name}_{num}_raw.tsv"
    filename_agg = f"cat2plot_{experiment}_{cat1name}_{cat2name}_{num}_agg.tsv"
    filename_plot = f"cat2plot_{experiment}_{cat1name}_{cat2name}_{num}"

    st.title(f'"{num}" vs "{cat1name}" and "{cat2name}"')

    data = pd.DataFrame(dict(
        num=numdata,
        cat1=cat1data,
        cat2=cat2data)).sort_values(by=['cat1', 'cat2'])

    agg = data.groupby(['cat1', 'cat2'])['num']\
            .agg(Mean=np.mean,
                 Median=np.median,
                 Stdev=np.std,
                 StdError=sem,)\
            .reset_index()
    
    fig = px.bar(agg, x='cat1', y='Mean', color='cat2',
                 barmode='group', error_y='StdError')
    
    fig.update_layout(                
        xaxis_title=cat1name,
        yaxis_title=f"{xtype}/{num}",
        legend=dict(title=cat2name),
    )

    plotly_config['toImageButtonOptions']['filename'] = filename_plot
    
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

