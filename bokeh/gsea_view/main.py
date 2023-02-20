from __future__ import generator_stop

import logging
import string
from functools import partial

from bokeh.layouts import column, row
from bokeh.models import (
    CheckboxGroup,
    ColumnDataSource,
    DataTable,
    ScientificFormatter,
    TableColumn,
)
from bokeh.models.widgets import Div, Select, TextAreaInput, TextInput
from bokeh.plotting import curdoc

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

w_selection_div = Div(text="")
w_geneset_filter = TextInput(title="Geneset filter", value="")
w_gene_list = TextAreaInput(title='Geneset genes',
                            rows=6, value="")
w_leadgene_list = TextAreaInput(title="Leading edge genes",
                                rows=6, value="")
w_newline = CheckboxGroup(labels=["Newline between genes"], active=[])

# Dataset
dataset_options = [
    (k, "{short_title}, {short_author}, {datatype}".format(**v))
    for k, v in datasets.items()]

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[0][0],
                             visible=False, height=30, width=400)


gsea_column_options = util.list2options(
    expset.get_gsea_columns(w_dataset_id.value))

w_diffexp = create_widget("diffexp", Select, title='Diff. Expression',
                          options=gsea_column_options, default='-')


gsea_sort_options = [('nes', 'NES'),
                     ('fdr', 'FDR'),
                     ('absnes', 'abs(NES)')]
w_sortorder = create_widget("sortorder", Select, title='Sort order',
                          options=gsea_sort_options, default='fdr')

w_debug = Div(text="")


def get_data():
    dataset_id, dataset = get_dataset()
    data = expset.get_gsea_data(dataset_id,
                                defield=w_diffexp.value,
                                gsfilter=w_geneset_filter.value,
                                sort_on=w_sortorder.value)
    studies = (data[['study_hash', 'study_title',
                     'study_author', 'study_year']]
               .drop_duplicates())

    # for easy reference - create short name for every study
    def nn(r):
        return r['study_author'].split(',')[0].split()[-1] + \
            " " + r['study_title'][:40] + \
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
    return data


source = ColumnDataSource(get_data())


w_table = DataTable(
    source=source,
    sizing_mode='stretch_both',
    margin=10,
    autosize_mode='fit_columns',
    index_position=None,
    columns=[
        TableColumn(field='geneset_title', title='Title', width=200),
        TableColumn(field='de_column', title='DE result',),
        TableColumn(field='study_short', title='Study', width=120,),
        TableColumn(field='no_genes',  width=40,
                    title='Genes in dataset'),
        TableColumn(field='no_lead_genes',  width=40,
                    title='No genes in Leading edge'),
        TableColumn(field='nes', title='NES', width=40,
                    formatter=ScientificFormatter(precision=2)),
        TableColumn(field='fdr', title='FDR', width=40,
                    formatter=ScientificFormatter(precision=1)),
    ])


def update_view(attr, old, new):
    curdoc().hold()
    dataset_id, dataset = get_dataset()

    data = get_data()
    source.data = data

    w_div_title_author.text = \
        f"""
            <b>{dataset['title']}</b>
            <i>{dataset['author']}</i><br>
            <b>Organism</b>: {dataset['organism']}<br>
            <b>Dataset</b>: {dataset['datatype']}
        """
    curdoc().unhold()


update_view(None, None, None)

w_diffexp.on_change("value", update_view)
w_sortorder.on_change("value", update_view)


def on_row_select(attrname, old, new):
    selectionIndex = source.selected.indices[0]

    def _get(f):
        return source.data[f][selectionIndex]

    if _get('geneset_type') == "geneset":
        gstype = "geneset"
    else:
        gstype = _get('geneset_type') + \
            f" ({_get('direction')})"

    w_selection_div.text = f"""
        <hr>
        <b>Differential expression rank:</b><br>
        {_get('de_column') }<br>

        <b>Against geneset </b>:<br>
        {_get('geneset_title') }<br>
        Type: { gstype }
        <br>
        <b>From:</b><br>
        {_get('study_title') }
        <i>{_get('study_author') }</i>
        ({_get('study_year') })
        <br><br>
        NES: {_get('nes') :.3g}<br>
        FDR: {_get('fdr') :.3g}<br>
    """

    w_gene_list.value = " ".join(_get("genes"))
    w_leadgene_list.value = " ".join(_get("lead_genes"))


def update_geneviews(attr, old, new):
    if len(w_newline.active) == 1:
        w_gene_list.value = w_gene_list.value.replace(' ', '\n')
        w_leadgene_list.value = w_gene_list.value.replace(' ', '\n')
    else:
        w_gene_list.value = w_gene_list.value.replace('\n', ' ')
        w_leadgene_list.value = w_gene_list.value.replace('\n', ' ')


selrow = column([
    w_selection_div, w_leadgene_list,
    w_gene_list,
    w_newline,
], sizing_mode="scale_width")


source.selected.on_change('indices', on_row_select)
w_newline.on_change('active', update_geneviews)
w_geneset_filter.on_change('value', update_view)


# Build the document
menucol = column(
    [row([w_dataset_id],
         sizing_mode="scale_width", ),
     row([w_div_title_author],
         sizing_mode="scale_width", ),
     row([w_geneset_filter], sizing_mode="scale_width", ),
     row([w_diffexp],
         sizing_mode="scale_width", ),
     row([w_sortorder],
         sizing_mode="scale_width", ),
     row([selrow], sizing_mode="stretch_both",),
     ],
    sizing_mode='fixed',
    width=350)

plotcol = column(
    [w_table], sizing_mode='stretch_both', )

curdoc().add_root(
    row([menucol, plotcol],
        sizing_mode='stretch_both'))
