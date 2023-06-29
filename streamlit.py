
from functools import partial
import pkg_resources


import streamlit as st
import plotly.express as px


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
        "format": "png",
        "width": 800,
        "height": 640,
        "scale": 3,
    })

st.sidebar.image(
    pkg_resources.resource_string('termite', 'data/termite.png')
)


@st.cache_data
def df_to_tsv(df):
    return df.to_csv(sep="\t").encode('utf-8')


view_placeholder = st.sidebar.empty()

experiment_data = db.get_experiments()
experiment_list = list(experiment_data['experiment'].unique())
experiment = util.selectbox_mem(st.sidebar, 'Experiment', experiment_list)
find_exp_id = experiment_data.query(f"experiment=='{experiment}'")
if len(find_exp_id) == 0:
    st.warning(f"Can not find experiment {experiment}")
    st.stop()
exp_id = find_exp_id.iloc[0]['exp_id']

gene_placeholder = st.sidebar.empty()

def categ_overview():
    import plotly.express as px
    all_catnames = db.get_obs_cat_names(exp_id)
    catname = st.sidebar.selectbox(
        "Categorical", all_catnames)
    st.title(f"Categorical: {catname}")
    cat_value_counts = db.get_obs_category(exp_id, catname)
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
        sql = sql.format(experiment=experiment)

    try:
        res = db.raw_sql(sql, raw=True)        
    except:
        raise
    
    else:
        st.write(res)

        try:
            st.download_button(
                "Download data",
                df_to_tsv(res),
                file_name="raw_sql_result.tsv",
                mime="text/tsv")
        except:
            st.write("Failed to prep for download")

def to_be_implemented():
    st.header("To be implemented")

    
subapps = {
    "Welcome": welcome.welcome,
    "Experiment": experiment_view.experiment_view,
    "Gene": to_be_implemented,
    "Metadata": categ_overview,
    "Run raw sql": run_sql,
    "Plot": 'plot',
    }


sapp = util.selectbox_mem(
    context = st.sidebar,
    label = 'Show',
    options = list(subapps.keys())
)

plotapps = {
    "Categorical": partial(categorical_plot,
                           experiment=experiment,
                           exp_id=exp_id,
                           plotly_config=plotly_config),
    "Two categoricals": partial(gene_vs_two_categories,
                                exp_id=exp_id,
                                experiment=experiment,
                                plotly_config=plotly_config),
    "Numerical": partial(numerical_numerical,
                         experiment=experiment,
                         exp_id=exp_id,
                         plotly_config=plotly_config),
}

if sapp != 'Plot':
    subapps[sapp]()
else:
    plotapp = util.selectbox_mem(
        context = st.sidebar,
        label = 'Plot',
        options = list(plotapps.keys())
    )
    plotapps[plotapp]()
    

