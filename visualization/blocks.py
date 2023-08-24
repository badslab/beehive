
from urllib.parse import quote
from textwrap import dedent

import streamlit as st
import pandas as pd


from termite.vis import util
from termite.vis import data as vdata


def get_dataset() -> pd.Series:
    search_string = st.sidebar.text_input("Search string")
    experiment = util.SqlSelect(
        "Experiment", key='exp',
        sql= f"""SELECT DISTINCT experiment, author, title, experiment
                   FROM dataset_md
                  WHERE ( title ILIKE '%{search_string}%'
                     OR author ILIKE '%{search_string}%'
                     OR experiment ILIKE '%{search_string}%' ) 
                   LIMIT 100 """)
    dataset = util.SqlSelect(
        "Dataset", key='dataset',
        oformat="{layertype} ({layername})",
        sql=f"""SELECT dataset, layername, layertype
                  FROM dataset_md
                 WHERE experiment='{experiment}'
                 LIMIT 100 """)
    dataset_rec = util.SqlOneRecord(
        sql=f"""SELECT *
                  FROM dataset_md
                 WHERE dataset='{dataset}' """)
    title(dataset_rec)
    return dataset_rec


def title(dataset_rec: pd.Series,
          context = None):
    if context is None:
        context = st

    ds = dataset_rec.loc['dataset']
    exp =  dataset_rec.loc['experiment']
    ds_quote = quote(ds)
    exp_quote = quote(exp)
    short_title = util.short_title(dataset_rec['title'])
    year = dataset_rec.loc['year']
    ds_url = f"/?app=dataset&dataset={ds_quote}&exp={exp_quote}"

    context.markdown(dedent(f"""
         <a href="{ds_url}" target="_self">
           <b>{short_title}</b>
           <i>{dataset_rec['author']}</i>
           ({year}, {dataset_rec['dataset']})</a>"""),
        unsafe_allow_html=True)

def get_expr(name, dataset_id):
    gene_name = util.SqlSelect(
        name, 
        sql=f"""SELECT DISTINCT gene
                  FROM help_gene
                 WHERE dataset_id='{dataset_id}'
                  ORDER BY gene""")
    gene = vdata.get_expr(dataset_id, gene_name)
    return gene_name, gene

def get_categorical(exp_id, name="Categorical", key="cat"):
    cat_name = util.SqlSelect(
        name, key=key,
        sql=f"""SELECT name, name, original_name
                  FROM help_obs
                 WHERE exp_id='{exp_id}'
                   AND dtype='categorical' """)
    cat = vdata.get_obs_cat(exp_id, cat_name)
    return cat_name, cat


def get_numerical(exp_id, name="Numerical", key="num",
                  exclude: list = []):
    
    if len(exclude) > 0:
        join = "','".join(exclude)
        exclude_stmt = f"AND name NOT IN ('{join}')"
    else:
        exclude_stmt = ""
     
    num_name = util.SqlSelect(
        name, key=key,
        sql=f"""SELECT name, name, original_name
                  FROM help_obs
                 WHERE exp_id='{exp_id}'
                   AND dtype IN ('int', 'float') {exclude_stmt}
        """)
    return num_name, vdata.get_obs_num(exp_id, num_name)


# def get_expr(dataset_id, gene):
#     return util.execute_sql(
#         sql=f"""SELECT obs, value
#                  FROM expr
#                 WHERE dataset_id = {dataset_id}
#                   AND gene = '{gene}' """).set_index('obs')['value']


def get_expr2(dataset_id, gene1, gene2):
    return util.execute_sql(
        sql=f"""SELECT a.value as g1, b.value as g2
                  FROM expr as a
                  JOIN expr as b on (a.obs=b.obs)
                  WHERE a.dataset_id={dataset_id}
                    AND b.dataset_id={dataset_id}
                    AND a.gene = '{gene1}'
                    AND b.gene='{gene2}'
             """)



def get_obs_num(exp_id, colname):
    return util.execute_sql(
        sql=f"""SELECT obs, value
                  FROM obs_num
                 WHERE exp_id={exp_id}
                   AND name='{colname}' """).set_index('obs')['value']


def get_dimred_name(experiment_id):
    return util.SqlSelect(
        "Dimred",
        sql=f"""
            SELECT DISTINCT dimred_name
              FROM help_obs
             WHERE exp_id={experiment_id}
               AND dtype='dimred' """)
