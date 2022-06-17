"""Simple expression visualization."""

from functools import partial
import logging
from pprint import pprint

import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource
from bokeh.models import DataTable, TableColumn, ScientificFormatter
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Select, TextInput, Div, Button, AutocompleteInput
from bokeh.plotting import figure, curdoc

from beehive import config, util, expset

lg = logging.getLogger('GeneExp')
lg.setLevel(logging.DEBUG)
lg.info("startup")



curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Gene Expression'

create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets()

args = curdoc().session_context.request.arguments

# WIDGETS
def_dataset_id = util.getarg(args, 'dataset_id',
                             list(datasets.keys())[0])

w_div_title_author = Div(text="")

dataset_options = [(k, "{title}, {author}".format(**v))
                   for k, v in datasets.items()]

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=def_dataset_id)

w_gene = create_widget("gene", AutocompleteInput, completions=[], default='APOE')
w_facet = create_widget("facet", Select, options=[], title="Group by")
w_plottype = create_widget("plottype", Select, title="Show",
                           options=["boxplot", "mean/std"])

w_download = Button(label='Download', align='end')
w_download_filename = Div(text="", visible=False,
                          name="download_filename")

# To display text if the gene is not found
w_gene_not_found = Div(text="")


#
# Data handling & updating interface
#
def get_genes():
    """Get available genes for a dataset."""
    dataset_id = w_dataset_id.value
    genes = sorted(list(expset.get_genes(dataset_id)))
    return genes


def get_facets():
    """Get available facets for a dataset."""
    dataset_id = w_dataset_id.value
    pprint(datasets[dataset_id])
    facets = sorted(list(datasets[dataset_id]['meta'].keys()))
    return facets


#
# Change & Initialize interface
#
def update_facets():
    """Update interface for a specific dataset."""
    facets = get_facets()
    w_facet.options = facets
    if w_facet.value not in facets:
        #set
        w_facet.value = \
            [f for f in facets
             if not f.startswith('_')][0]

def update_genes():
    """Update genes widget for a dataset."""
    genes = get_genes()
    w_gene.completions = genes
    if w_gene.value not in genes:
        if 'APOE' in genes:
            w_gene.value = 'APOE'
        else:
            w_gene.value = genes[0]


update_facets()
update_genes()


def get_data() -> pd.DataFrame:
    """Retrieve data from a dataset, gene & facet."""
    dataset_id = w_dataset_id.value
    data = expset.get_gene_meta_agg(
        dsid=dataset_id,
        gene=w_gene.value,
        meta=w_facet.value)

    data['perc'] = 100 * data['count'] / data['count'].sum()

    # default settings for a boxplot -
    # override if other further down

    data['_segment_top'] = data['q99']
    data['_bar_top'] = data['q75']
    data['_bar_median'] = data['median']
    data['_bar_bottom'] = data['q25']
    data['_segment_bottom'] = data['q01']

    return data


def get_dataset() -> dict:
    """Return the current dataset record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]


#
# Create plot
#
plot = figure(background_fill_color="#efefef", x_range=[],
              plot_height=400, title="Plot",
              toolbar_location='right')
source = ColumnDataSource(get_data())
table = DataTable(source=source,
                  margin=10,
                  index_position=None,
                  columns=[
                      TableColumn(field='cat_value', title='Category'),
                      TableColumn(field='count', title='Count',
                                  formatter=ScientificFormatter(precision=0)),
                      TableColumn(field='perc', title='Percentage',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='mean', title='Mean',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='median', title='Median',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q01', title='1% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q25', title='20% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q75', title='80% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q99', title='99% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                  ])


# create segments
elements = dict(
    vbar=plot.vbar(source=source, x='cat_value', top='_bar_top',
                   bottom='_bar_bottom', width=0.85, name="barplot",
                   fill_color=config.colors.color6,
                   line_color="black"),
    seg_v_up=plot.segment(source=source, x0='cat_value', x1='cat_value',
                          y0='_bar_top', y1='_segment_top', line_color='black'),
    seg_h_up=plot.rect(source=source, x='cat_value', height=0.001,
                       y='_segment_top', width=0.4, line_color='black'),
    seg_v_dw=plot.segment(source=source, x0='cat_value', x1='cat_value',
                          y0='_segment_bottom', y1='_bar_bottom',
                          line_color='black'),
    seg_h_dw=plot.rect(source=source, x='cat_value', height=0.001,
                       y='_segment_bottom', width=0.4, line_color='black'),
    seg_h_med=plot.rect(source=source, x='cat_value', height=0.001,
                        y='_bar_median', width=0.85, line_width=2,
                        line_color='black'),
)


def update_plot():
    """Update the plot."""

    global plot, source
    data = get_data()
    dataset_id, dataset = get_dataset()
    facet = w_facet.value
    gene = w_gene.value

    w_div_title_author.text = \
        f"""
        <ul>
          <li><b>Title:</b> {dataset['title']}</li>
          <li><b>Author:</b> {dataset['author']}</li>
        </ul>
        """

    gene = w_gene.value
    plot.x_range.factors = list(sorted(data['cat_value']))
    w_download_filename.text = f"exp_{dataset_id}_{facet}_{gene}.tsv"

    if w_plottype.value == "boxplot":
        pttext = 'Boxplot'
    else:
        data['_segment_top'] = data['mean'] + data['std']
        data['_bar_top'] = data['mean']
        data['_bar_bottom'] = 0
        data['_segment_bottom'] = 0
        data['_bar_median'] = data['mean']
        pttext = 'Mean/std'


    ymax = max(1, 1.01 * data['_segment_top'].max())
    source.data = data

    title = dataset['title'][:60]
    plot.title.text = (f"{pttext} {gene} vs {facet} - "
                       f"({dataset_id}) {title}...")
    print("setting ymax to", str(ymax)[:50])
    plot.y_range.update(end = ymax, start=-0.1)


update_plot()


#
# widget callbacks
#
def _dataset_change(attr, old, new):
    """Dataaset change."""
    curdoc().hold()
    update_facets()
    update_genes()
    update_plot()
    curdoc().unhold()


def _update_plot(attr, old, new):
    update_plot()


download_callback = CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")


w_gene.on_change("value", _update_plot)
w_dataset_id.on_change("value", _dataset_change)
w_facet.on_change("value", _update_plot)
w_plottype.on_change("value", _update_plot)
w_download.js_on_click(download_callback)


#
# Build the document
#
curdoc().add_root(
    column([
        row([w_dataset_id, w_download_filename],
            sizing_mode='stretch_width'),
        row([w_gene, w_facet, w_plottype, w_download],
            sizing_mode='stretch_width'),
        row([w_div_title_author],
            sizing_mode='stretch_width'),
        row([w_gene_not_found],
            sizing_mode='stretch_width'),
        row([plot],
            sizing_mode='stretch_width'),
        row([table],
            sizing_mode='stretch_width'),
    ], sizing_mode='stretch_width')
)
