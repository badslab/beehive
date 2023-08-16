
from textwrap import dedent
from urllib.parse import quote

import pandas as pd
import streamlit as st

import views

from termite.vis import data as vdata

import blocks


def dataset():

    dataset_rec = blocks.get_dataset()  
    ds_id =  dataset_rec.loc['dataset_id']
    ds = dataset_rec.loc['dataset']
    exp_id =  dataset_rec.loc['experiment_id']
    exp =  dataset_rec.loc['experiment']
    ds_quote = quote(ds)
    exp_quote = quote(exp)
    
    ds_stats = vdata.get_dataset_stats(ds_id)
    del ds_stats['dataset_id']

    md = vdata.get_obs_data_exp(exp_id)
    
    st.title(dataset_rec['title'])
    st.markdown(dedent(
        """
        *{author}*, {year}
        
        **dataset**: {dataset} ({dataset_id})
        from **experiment**: {experiment} ({experiment_id})
        and **study** {study}
        
        * **organism**: {organism}
        * **Layer**: {layername} ({layertype})
        * **Total datapoints**: {no_datapoints:_g}
        * **No Observations (cells)**: {no_cells:_g}
        * **No Variables (genes)**: {no_genes:_g}
        """.format(**dataset_rec, **ds_stats)))

    topgenes = vdata.get_topgenes(ds_id)
    views.BarPlot(data=topgenes[['gene', 'sumval']],
                  subject='Most higly expressed genes').view()

    def format_gene(gene):
        gene_q = quote(gene)
        url = f"/?app=gene&dataset={ds_quote}&exp={exp_quote}&gene={gene_q}"
        return f'<a href="{url}" target="_self">{gene}</a>'
    def format_cat(c):
        cq = quote(c)
        url = f"/?app=categ&dataset={ds_quote}&exp={exp_quote}&categorical={cq}"
        return f'<a href="{url}" target="_self">{c}</a>'
    def format_num(c):
        cq = quote(c)
        url = f"/?app=numer&dataset={ds_quote}&exp={exp_quote}&numerical={cq}"
        return f'<a href="{url}" target="_self">{c}</a>'
    
    obscat = "&nbsp;&middot;&nbsp;".join(
        map(format_cat, md[md['dtype'] == "categorical"]['name']))
    obsint = "&nbsp;&middot;&nbsp;".join(
        map(format_num, md[md['dtype'] == "int"]['name']))
    obsfloat = "&nbsp;&middot;&nbsp;".join(
        map(format_num, md[md['dtype'] == "float"]['name']))
    embed_gene = "&nbsp;&middot;&nbsp;".join(
        map(format_gene, topgenes['gene']))
    obsdim = '`' + "`, `".join(md[md['dtype'] == "dimred"]['name']) + '`'

    st.markdown(dedent(
        """
        * **Top Genes**: {embed_gene}
        * **Categorical variables**: {obscat}
        * **Integer variables**: {obsint}
        * **Float variables**: {obsfloat}
        * **Dim.red variables**: {obsdim}
        """.format(obscat=obscat, obsdim=obsdim, embed_gene=embed_gene,
                   obsfloat=obsfloat, obsint=obsint)),
                unsafe_allow_html=True)


def categorical():
    dataset_rec = blocks.get_dataset()
    categ_name = blocks.get_categorical(dataset_rec['experiment_id'])
    cat = vdata.get_obs_cat(dataset_rec['experiment_id'], categ_name)
    cvc = cat.value_counts()
    cvc.index = cvc.index.astype(str)
    cvc.index.name = categ_name
    cvc = pd.DataFrame(cvc).reset_index()
    views.BarPlot(data=cvc, subject=categ_name).view()
    
def catcat():
    dataset_rec = blocks.get_dataset()
    cat_name_1 = blocks.get_categorical(dataset_rec['experiment_id'],
                                        name='Categorical 1', key='cat')
    cat_name_2 = blocks.get_categorical(dataset_rec['experiment_id'],
                                        name='Categorical 2', key='cat2')
    cat1 = vdata.get_obs_cat(dataset_rec['experiment_id'], cat_name_1)
    cat2 = vdata.get_obs_cat(dataset_rec['experiment_id'], cat_name_2)

    
    data = pd.DataFrame({cat_name_1:cat1, cat_name_2:cat2})\
        .pivot_table(columns=cat_name_1, index=cat_name_2, aggfunc=len)
    
    views.HeatMap(data=data, subject=f"{cat_name_1} vs {cat_name_2}").view()
    

    
def numerical():
    dataset_rec = blocks.get_dataset()
    num_name = blocks.get_numerical(dataset_rec['experiment_id'])
    num = vdata.get_obs_num(dataset_rec['experiment_id'], num_name)
    vargs = dict(data=num, subject=num_name)
    views.Multi(
        views.Ecdf(**vargs),
        views.Histogram(**vargs)
        ).view()
