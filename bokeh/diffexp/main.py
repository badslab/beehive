"""Simple diff. expression visualization."""

from functools import partial
import logging
import math
from pprint import pprint

import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource
from bokeh.models import DataTable, TableColumn, ScientificFormatter
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Select, TextInput, Div, Button
from bokeh.plotting import figure, curdoc
from bokeh.models import HTMLTemplateFormatter

import colorcet as cc

from beehive import config, util, expset


lg = logging.getLogger('DiffExp')
lg.setLevel(logging.DEBUG)
lg.info("startup")


curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Differential Expression'

create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets(has_de=True)

lg.info(f"Discovered {len(datasets)} datasets")

args = curdoc().session_context.request.arguments

# WIDGETS
w_div_title_author = Div(text="")


def_dataset_id = util.getarg(args, 'dataset_id',
                             list(datasets.keys())[0])
print(def_dataset_id)


dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=def_dataset_id,
                             visible=False,)

w_gene = create_widget("gene", TextInput, title="Genes (space separated)",
                       default='APOE TREM2 PLP1 GFAP RBFOX3')

w_facet = create_widget("facet", Select, options=[],
                        title="Differential expression across:")

siblings = expset.get_dataset_siblings(w_dataset_id.value)

sibling_options = []
for k, v in siblings.items():
    sname = f"{v['organism']} / {v['datatype']}"
    sibling_options.append((k, sname))
    print(k, sname, def_dataset_id)

w_sibling = create_widget("view", Select,
                          options=sibling_options,
                          default=def_dataset_id,
                          update_url=False)


#w_download = Button(label='Download', align='end')
w_problem = Div(text="", visible=True)


#
# Data handling & updating interface
#

def get_dataset():
    """Return the current dataset record."""
    dataset_id = w_dataset_id.value.strip()
    if not dataset_id:
        dataset_id = def_dataset_id
    lg.info("Found dataset: " + dataset_id)
    return dataset_id, datasets[dataset_id]


def get_genes():
    """Get available genes for a dataset."""
    dataset_id, dataset = get_dataset()
    lg.info("######" + dataset_id)
    genes = sorted(list(expset.get_genes(dataset_id)))
    return genes


def get_facets():
    """Get available facets for a dataset."""
    dataset_id, dataset = get_dataset()
    lg.info(f"Get facets for {dataset_id}")
    facets = expset.get_defields(dataset_id)
    facets = sorted(list(facets))
    lg.info(f"Found facets {facets}")
    return facets

#
# Change & Initialize interface
#


def update_facets():
    """Update interface for a specific dataset."""

    facets = get_facets()

    w_facet.options = facets
    if w_facet.value not in facets:
        w_facet.value = \
            [f for f in facets
             if not f.startswith('_')][0]
        print('-' * 30 + str(list(facets)))


# Ensure the facet widget is up to date
update_facets()


def get_current_genes():
    """Get relevant set of genes."""
    dataset_id, dataset = get_dataset()
    genes = w_gene.value
    genes = genes.split()

    allgenes = set(expset.get_genes(dataset_id))
    genes = set(genes) & allgenes
    genes = list(sorted(genes))

    if len(genes) == 0:
        dsid, ds = get_dataset()
        topgenes = ds['diffexp'][w_facet.value]['topgenes']
        w_problem.text = "No genes found - setting top genes"
        w_gene.value = " ".join(sorted(topgenes))

    return genes


def get_data():
    """Retrieve DE data for a dataset, gene & facet."""

    dataset_id = w_dataset_id.value
    genes = get_current_genes()

    lg.info(f"Get DE data for {dataset_id} / {w_facet.value}")
    if len(genes) == 0:
        return None, None

    data = expset.get_dedata(dataset_id, w_facet.value, genes)

    print(data)

    lfc = data.xs('lfc', axis=1, level=1)
    padj = data.xs('padj', axis=1, level=1)
    padj.columns = padj.columns + "__padj"

    notsig = lfc > 0.05
    data = pd.concat([lfc, padj], axis=1)

    try:
        # see if the index is numerical -
        # if so sort along that order
        data.sort_index(level=0, key=lambda x: x.astype(float))
    except:
        data = data.sort_index()

    vmax = (lfc * notsig).abs().max().max()
    data = data.reset_index()
    return vmax, data


#
# Create plot
#
source = ColumnDataSource(dict(cat_value=[]))
table = DataTable(source=source,
                  margin=10,
                  index_position=None,
                  columns=[
                      TableColumn(field='cat_value', title='Category'),
                  ])


def get_html_formatter(col):
    template = f"""<%= {col}__html %>"""
    return HTMLTemplateFormatter(template=template)


def update_table():
    lg.info("Updating table")

    vmax, data = get_data()
    genes = get_current_genes()
    lg.info("UT: " + " ".join(genes))

    ds = get_dataset()

    if data is None:
        w_problem.text = "No genes found?"
        return

    def formatHtml(row, gene):
        padj = row[f"{gene}__padj"]
        lfc = row[f"{gene}"]
        lfc = max(-vmax, min(vmax, lfc))

        lp = 99 if padj == 0 else int(min(99, -math.log10(padj)))

        lfc_col_no = int(255 * (lfc + vmax) / (2 * vmax))
        lfc_col = cc.CET_D3[lfc_col_no]
        lfc_col_txt = "#000"

        p_col_no = min(255, int(255 * (lp / 99)))
        p_col = cc.blues[p_col_no]
        p_col_txt = '#000'

        if lp > 75:
            p_col_txt = '#FFF'

        if padj > 0.05:
            lfc_col = config.colors.lightgrey
            lfc_col_txt = config.colors.darkgrey
            p_col = config.colors.lightgrey
            p_col_txt = config.colors.lightgrey

        return f"""
           <span class="lfc" style="color:{lfc_col_txt};background-color:{lfc_col};">
               {lfc:.1f}</span>
           <span class="pval" style="color:{p_col_txt};background-color:{p_col};">
               {lp}</span>
        """.strip()

    for g in genes:
        data[f"{g}__html"] = data.apply(formatHtml, axis=1, gene=g)

    newcols = [TableColumn(field='cat_value', title='Category')]

    for gene in genes:
        newcols.append(
            TableColumn(field=gene, title=gene,
                        formatter=get_html_formatter(gene)))

    # update the doc in one go...
    curdoc().hold()
    table.columns = newcols
    source.data = data
    curdoc().unhold()


update_table()

#
# widget callbacks
#


def cb_dataset_change(attr, old, new):
    """Dataaset change."""
    lg.info("CB dataset change")
    curdoc().hold()
    update_facets()
    update_table()
    curdoc().unhold()


def cb_genes_change(attr, old, new):
    """Gene change."""
    update_table()


def cb_facet_change(attr, old, new):
    lg.info("CB facet change")
    update_table()


w_gene.on_change("value", cb_genes_change)
w_dataset_id.on_change("value", cb_dataset_change)
w_facet.on_change("value", cb_facet_change)


#
# Build the document
#
curdoc().add_root(
    column([
        row([w_dataset_id, w_sibling, w_facet],
            sizing_mode='stretch_width'),
        row([w_gene],
            sizing_mode='stretch_width'),
        row([w_div_title_author],
            sizing_mode='stretch_width'),
        row([w_problem],
            sizing_mode='stretch_width'),
        row([table],
            sizing_mode='stretch_width'),
    ], sizing_mode='stretch_width')
)
