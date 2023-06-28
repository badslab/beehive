
from functools import partial
import pkg_resources


import streamlit as st
import plotly.express as px


from streamlit_option_menu import option_menu


from termite import db
from termite.streamlit import welcome, util, experiment_view
from termite.streamlit.categorical_plot import categorical_plot
from termite.streamlit.gene_vs_two_categories import gene_vs_two_categories
from termite.streamlit.numerical_numerical import numerical_numerical


st.set_page_config(layout="wide", page_title="Termite")
px.defaults.color_discrete_sequence = px.colors.qualitative.T10

plotly_config = dict(
    displaylogo = False,
    toImageButtonOptions = {
        "format": "svg",
    })

st.sidebar.image(
    pkg_resources.resource_string('termite', 'data/termite.png')
)


@st.cache_data
def df_to_tsv(df):
    return df.to_csv(sep="\t").encode('utf-8')


view_placeholder = st.sidebar.empty()

experiment_data = db.get_exp_datatypes()
experiment_list = list(experiment_data['experiment'].unique())

experiment = util.selectbox_mem(st.sidebar, 'Experiment', experiment_list)

datatypes = list(experiment_data\
    [experiment_data['experiment'] == experiment]\
    ['datatype'].unique())

datatype = util.selectbox_mem(st.sidebar, 'Datatype', datatypes)


gene_placeholder = st.sidebar.empty()


def categ_overview():
    import plotly.express as px
    all_catnames = db.get_obs_cat_names(experiment)
    catname = st.sidebar.selectbox(
        "Categorical", all_catnames)
    st.title(f"Categorical: {catname}")
    cat_value_counts = db.get_obs_category(experiment, catname)
    fig = px.bar(cat_value_counts, x='value', y='count')
    st.plotly_chart(fig, config=plotly_config)
    
    st.download_button(
        "Download data",
        df_to_tsv(cat_value_counts),
        file_name=f"category_{experiment}_{catname}.tsv",
        mime="text/tsv")


def run_sql():
    st.title("Raw SQL")
    sql = st.text_area("SQL:", "SELECT *\nFROM expr\nLIMIT 5")
    if '{' in sql:
        sql = sql.format(experiment=experiment, datatype=datatype)
    try:
        res = db.raw_sql(sql)
    except:
        raise
    else:
        st.write(res)
        st.download_button(
            "Download data",
            df_to_tsv(res),
            file_name="raw_sql_result.tsv",
            mime="text/tsv")


subapps = {
    "Welcome": welcome.welcome,
    "Experiment": experiment_view.experiment_view,
    "Metadata": categ_overview,
    "Categorical Plot": partial(
        categorical_plot,
        experiment=experiment,
        datatype=datatype,
        plotly_config=plotly_config),
    "Gene / two categoricals": partial(gene_vs_two_categories,
                                       experiment=experiment,
                                       datatype=datatype,
                                       plotly_config=plotly_config),
    "Numerical": partial(numerical_numerical,
                         experiment=experiment,
                         datatype=datatype,
                         plotly_config=plotly_config),
    "Run Sql": run_sql,
}

# Horizontal menu
# see - https://github.com/victoryhb/streamlit-option-menu
with view_placeholder:
    sapp = option_menu(
        None, list(subapps.keys()),
        menu_icon="cast",
        default_index=0,
        orientation="vertical")
        
subapps[sapp]()

