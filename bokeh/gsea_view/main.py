from __future__ import generator_stop

import logging
import string
from functools import partial
from os.path import dirname, join

import pandas as pd
from bokeh.layouts import column, row
from bokeh.models import (
    ColumnDataSource,
    DataTable,
    NumberFormatter,
    Range1d,
    ScientificFormatter,
    TableColumn,
)
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import AutocompleteInput, Button, Div, Select
from bokeh.plotting import curdoc, figure
from bokeh.transform import CategoricalColorMapper, jitter

import beehive.exceptions as bex
from beehive import config, expset, util

lg = logging.getLogger('GseaView')
lg.setLevel(logging.DEBUG)
lg.info("startup")

VIEW_NAME = "gsea_table"

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'GSEA Results'

create_widget = partial(util.create_widget, curdoc=curdoc())

args = curdoc().session_context.request.arguments

datasets = expset.get_datasets(view_name=VIEW_NAME)


def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]


# @lru_cache(16)
def get_diffexp_field_options():
    dataset_id = w_dataset_id.value
    defields = expset.get_defields(dataset_id)
    defields = [(f"{a}__lfc", b) for (a,b) in defields]
    return defields


# WIDGETS
w_div_title_author = Div(text="")
warning_div = Div(text="nothing good about this",
                  sizing_mode="stretch_both",
                  css_classes=['beehive_warning'])

# Dataset
dataset_options = [
    (k, "{short_title}, {short_author}, {datatype}".format(**v))
    for k, v in datasets.items()]

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[0][0],
                             visible=False, height=30, width=400)

_dexo = get_diffexp_field_options()
w_diffexp = create_widget("diffexp", Select, title='Diff. Expression',
                          options=_dexo, default=_dexo[0][0])

w_debug = Div(text="")


def get_data():
    dataset_id, dataset = get_dataset()
    defield = w_diffexp.value
    data = expset.get_gsea_data(dataset_id, defield)
    print(data.head(2).T)
    studies = (data[['study_hash', 'study_author', 'study_year']]
               .drop_duplicates())
    print(studies)

    # for easy reference - create short name for every study
    def nn(r):
        return r['study_author'].split(',')[0].split()[-1] + \
            f" {r['study_year']}"
    studies['study_short'] = studies.apply(nn, axis=1)
    i = -1

    def get_letter():
        nonlocal i
        i += 1
        return string.ascii_lowercase[i]

    studies['study_short'][studies['study_short'].duplicated()] \
        = studies['study_short'][studies['study_short'].duplicated()].apply(
            lambda x: x + " " + get_letter())
    data = data.merge(studies[['study_hash', 'study_short']], how='left',
                      on='study_hash')
    print(data.head(2).T)
    return data

source = ColumnDataSource(get_data())


w_table = DataTable(
    source=source,
    sizing_mode='stretch_both',
    margin=10,
    index_position=None,
    columns=[
        TableColumn(field='title', title='Title'),
        TableColumn(field='study_short', title='Study'),
        TableColumn(field='no_genes', title='Genes in dataset'),
        TableColumn(field='nes', title='NES',
                    formatter=ScientificFormatter(precision=2)),
        TableColumn(field='fdr', title='FDR',
                    formatter=ScientificFormatter(precision=1)),
    ])


def update_view(attr, old, new):
    print('update')
    curdoc().hold()
    dataset_id, dataset = get_dataset()
    defield = w_diffexp.value

    w_debug.text = str(defield)
    lg.info(f"DE Field {defield}")

    data = get_data()
    source.data = data

    w_div_title_author.text = \
        f"""
        <ul>
          <li><b>Title:</b> {dataset['title']}</li>
          <li><b>Author:</b> {dataset['author']}</li>
          <li><b>Organism / Datatype:</b>
              {dataset['organism']} / {dataset['datatype']}</li>
        </ul>
        """
    print('done')
    curdoc().unhold()


update_view(None, None, None)

w_diffexp.on_change("value", update_view)


# Build the document
menucol = column(
    [row([w_dataset_id],
         sizing_mode="scale_width", ),
     row([w_div_title_author],
         sizing_mode="scale_width", ),
     row([w_diffexp],
         sizing_mode="scale_width", ),
    ],
    sizing_mode='fixed',
    width=350)

plotcol = column(
    [w_table], sizing_mode='stretch_both', )


curdoc().add_root(
    row([menucol, plotcol],
    sizing_mode='stretch_both'))
