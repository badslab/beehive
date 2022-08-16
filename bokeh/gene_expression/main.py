"""Simple boxplot gene expression visualization."""

# Notes:
#
# I removed the plottype for the time being - mean/std plots
# are strange when the min value is not zero - some sets
# have negative (log) values - discuss if we can fix this
# - or simply do not show mean/std plots? - solution for now.
#


from functools import partial
import logging
from pprint import pprint

import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models import DataTable, TableColumn, ScientificFormatter
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, TextInput, Div,
                                  Button, AutocompleteInput)
from bokeh.plotting import figure, curdoc

from beehive import config, util, expset

from bokeh.transform import jitter

lg = logging.getLogger('GeneExp')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Gene/Protein Expression'

create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets()

args = curdoc().session_context.request.arguments

# WIDGETS
w_div_title_author = Div(text="")

# Dataset

dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[0][0],
                             visible=True,height = 10, width = 100)

w_sibling = create_widget("view", Select,
                          options=[],
                          default=w_dataset_id.value,
                          update_url=False,height = 20, width = 150)


def update_sibling_options():
    siblings = expset.get_dataset_siblings(w_dataset_id.value)
    sibling_options = []
    # check all organisms
    organisms = []
    for k, v in siblings.items():
        organisms = organisms + [v['organism']]

    if len(list(set(organisms))) == 1:
        for k, v in siblings.items():
            # sname = f"{v['organism']} / {v['datatype']}"
            sname = f"{v['datatype']}"
            sibling_options.append((k, sname))
    else:
        for k, v in siblings.items():
            sname = f"{v['organism']} / {v['datatype']}"
            sibling_options.append((k, sname))

    w_sibling.options = sibling_options


update_sibling_options()


w_gene = create_widget("gene", AutocompleteInput,
                       completions=[], default='APOE', case_sensitive=False, height = 50, width = 150)
w_facet = create_widget("facet", Select, options=[], title="Group by level 1:",height = 50, width = 150)
w_facet2 = create_widget("facet2", Select, options=[], title="Group by level 2:",height = 50, width = 150)
# w_facet3 = create_widget("facet3", Select, options=[], title="Show class of avg expression in:",height = 50, width = 150)
w_download = Button(label='Download', align='end', height = 50, width = 150)
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


def update_facets():
    """Update interface for a specific dataset."""
    options = expset.get_facet_options(w_dataset_id.value,only_categorical = True)
    # options_with_skip = expset.get_facet_options(w_dataset_id.value,only_categorical = True, include_skip = True)
    w_facet.options = options
    w_facet2.options = options
    # w_facet3.options = options_with_skip

    if w_facet.value not in [x[0] for x in options]:
        # set a default
        w_facet.value = options[0][0]

    w_facet2.options = list(filter(lambda x: x[0] != w_facet.value,w_facet2.options))
    if w_facet2.value not in [x[0] for x in options]:
        # set a default
        w_facet2.value = options[0][0]
    
    w_facet.options = list(filter(lambda x: x[0] != w_facet2.value,w_facet.options))

    # if w_facet3.value not in [x[0] for x in options_with_skip]:
    #     # set a default
    #     w_facet3.value = options_with_skip[0][0]


def facet_groups(facet):
    result_dict = expset.get_facet_groups(w_dataset_id.value, facet)
    return result_dict


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
    gene = w_gene.value
    facet = w_facet.value
    facet2 = w_facet2.value
    # facet3 = w_facet3.value

    # TODO get the order from here and sort the data later..
    # will also get color
    # facet_groups = facet_groups(w_facet.value)

    # ==> check if there is full order
    # if yes, call get_gene_meta_agg with order = True
    # can order with the list of orders in expset.get_gene_meta_agg before function returns

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene}")


    data = expset.get_gene_meta_three_facets(dataset_id,gene,facet,facet2,"mouse.id")
    data_no_dups = data.drop_duplicates("cat_value")

    return data,data_no_dups


def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]


#
# Create plot
#
plot = figure(background_fill_color="#efefef", x_range=[],title="Plot",
              toolbar_location='right', tools="save", sizing_mode = "scale_both")

data,data_no_dups = get_data()
source = ColumnDataSource(data)
source_no_dups = ColumnDataSource(data_no_dups)
table = DataTable(source=source_no_dups,
                  margin=10,
                  index_position=None,
                  sizing_mode = "scale_both",
                  columns=[
                      TableColumn(field='cat_value', title='Category'),
                      TableColumn(field='count', title='No Samples/Cells',
                                  formatter=ScientificFormatter(precision=0)),
                      TableColumn(field='perc', title='% Samples/Cells',
                                  formatter=ScientificFormatter(precision=1)),
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

# meta3 = w_facet3.value
meta3 = "mouse.id"
# create plot elements - these are the same for boxplots as mean/std type plots
elements = dict(
    vbar=plot.vbar(x="cat_value", top='_bar_top',
                bottom='_bar_bottom', source = source, width=0.85, name="barplot",
                fill_color=config.colors.color6,
                line_color="black"),
    seg_v_up=plot.segment(source=source, x0='cat_value', x1='cat_value',
                          y0='_bar_top', y1='_segment_top',
                          line_color='black'),
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
    jitter_points = plot.scatter(x=jitter('cat_value', width=0.4, range=plot.x_range), y=f'mean_{meta3}', size=5, alpha=0.4, source=source,legend_label = f"{meta3}")

)

X_AXIS_LABELS_ORIENTATION = 3.14/2

yspacer = (data['_segment_top'].max() - data['_segment_bottom'].min()) / 20

ymax = data['_segment_top'].max() + yspacer
ymin = data['_segment_bottom'].min() - yspacer

plot.update(y_range=Range1d(ymin, ymax))
plot.xaxis.major_label_orientation = X_AXIS_LABELS_ORIENTATION


def cb_update_plot(attr, old, new):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, source

    data,data_no_dups = get_data()
    update_facets()
    dataset_id, dataset = get_dataset()
    facet = w_facet.value
    facet2 = w_facet2.value

    gene = w_gene.value
    # meta3 = w_facet3.value
    w_div_title_author.text = \
        f"""
        <ul>
          <li><b>Title:</b> {dataset['title']}</li>
          <li><b>Author:</b> {dataset['author']}</li>
          <li><b>Organism / Datatype:</b>
              {dataset['organism']} / {dataset['datatype']}</li>
        </ul>
        """

    gene = w_gene.value
    plot.x_range.factors = sorted(list(set(data["cat_value"])),key=lambda tup: tup[0])
    w_download_filename.text = f"exp_{dataset_id}_{facet}_{gene}.tsv"

    source.data = data
    source_no_dups.data = data_no_dups
    # plan for 5% space above & below
    yspacer = (data['_segment_top'].max() - data['_segment_bottom'].min()) / 20

    ymax = data['_segment_top'].max() + yspacer
    ymin = data['_segment_bottom'].min() - yspacer

    lg.warning(f"## Y MIN/MAX {ymin:.2g} / {ymax:.2g} ")
    plot.y_range.update(start=ymin, end=ymax)

    title = dataset['short_title']
    if len(title) > 80:
        title = title[:77] + '...'
    plot.title.text = (f"Boxplot of {gene} vs {facet}"
                       f" - {dataset['organism']}"
                       f" - {dataset['first_author']} - {title}")
    plot.yaxis.axis_label = f"{dataset['datatype']}"
    plot.xaxis.major_label_orientation = X_AXIS_LABELS_ORIENTATION

    curdoc().unhold()


# convenience shortcut
update_plot = partial(cb_update_plot, attr=None, old=None, new=None)

# run it directly to ensure there are initial values
update_plot()



def cb_dataset_change(attr, old, new):
    """Dataset change."""
    lg.info("Dataset Change to:" + new)
    curdoc().hold()
    update_facets()
    update_genes()
    update_sibling_options()
    update_plot()


def cb_sibling_change(attr, old, new):
    lg.debug("Sibling change: " + new)
    w_dataset_id.value = new
    update_plot()


cb_download = CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")


w_gene.on_change("value", cb_update_plot)
w_sibling.on_change("value", cb_sibling_change)
w_dataset_id.on_change("value", cb_dataset_change)
w_facet.on_change("value", cb_update_plot)
w_download.js_on_click(cb_download)
w_dataset_id.on_change("value", cb_update_plot)
w_facet2.on_change("value",cb_update_plot)
# w_facet3.on_change("value",cb_update_plot)

#
# Build the document
#
curdoc().add_root(row([
    column([
        column([
        row([w_gene, w_facet,w_facet2],sizing_mode='scale_both'),
        # row([w_sibling, w_download,w_facet3],sizing_mode='scale_both'),
        row([w_sibling, w_download],sizing_mode='scale_both'),

        ]),
        column([w_div_title_author], sizing_mode='scale_both'),
        column([w_dataset_id],sizing_mode='scale_both'),
        column([table],sizing_mode='scale_both')
        ]),
    column([
        column([plot], sizing_mode='scale_both')
    ])
], sizing_mode='scale_both')
)
