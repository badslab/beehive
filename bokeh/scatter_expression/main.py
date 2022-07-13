import logging
from functools import partial
import logging
from pprint import pprint
import numpy as np
import pandas as pd
import polars as pl

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource
from bokeh.models import DataTable, TableColumn, ScientificFormatter
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, TextInput, Div,
                                  Button, AutocompleteInput)
from bokeh.plotting import figure, curdoc, show
from bokeh.transform import factor_cmap, factor_mark
from bokeh.palettes import Category10

from beehive import config, util, expset

lg = logging.getLogger('ScatterExpression')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Scatter Expression'


create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets()

args = curdoc().session_context.request.arguments



w_div_title_author = Div(text="")

#default dataset... ex: m.how1m.2
def_dataset_id = util.getarg(args, 'dataset_id',
                             list(datasets.keys())[5])
######SETTING IT MANUALLY FOR NOW############
def_dataset_id = "h.man2m.1"

#datasets with titles
dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=def_dataset_id,
                             visible=False,)

#keys are the dataset names, and values are
#'author', 
# 'organism',  mouse/human
# 'datatype', ex: protein concentration
# 'description', 
# 'diffexp' : (mouse_models': {'keys': ['BV2'], 'topgenes': ['H2AFY', 'LMNB2', 'H1F0', 'CST3', 'DPYSL2', 'CKB', 'MECP2', 'LDHB'])
#  'dimred' ??
# 'meta', {'model': {'dtype': 'categorical'}} <= to be looked at
#  'short_title', 
# 'slug', 'title', 'year', 'study': how1m..how2m...
#  'short_author'])

## note: human/plant has meta as 'models' , mouse has meta as 'model'
siblings = expset.get_dataset_siblings(w_dataset_id.value)


#similar, create the options of the siblings datasets
sibling_options = []
for k, v in siblings.items():
    sname = f"{v['organism']} / {v['datatype']}"
    sibling_options.append((k, sname))
    print(k, sname, def_dataset_id)


#make a widget for the siblings
w_sibling = create_widget("view", Select,
                          options=sibling_options,
                          default=def_dataset_id,
                          update_url=False)

#making widgets for the following:

#gene, has an autocomplete feature.. write in text
w_gene1 = create_widget("gene1", AutocompleteInput,
                       completions=[], default='APOE')

w_gene2 = create_widget("gene2", AutocompleteInput,
                       completions=[], default='APOE')

#drop down menu to group by items
w_facet = create_widget("facet", Select, options=[], title="Group by")

#download button..
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

##only model ??
def get_facets():
    """Get available facets for a dataset."""
    dataset_id = w_dataset_id.value
    facets = sorted(list(datasets[dataset_id]['meta'].keys()))
    return facets

#
# Change & Initialize interface
#


def update_genes():
    """Update genes widget for a dataset."""
    genes = get_genes()
    w_gene1.completions = genes
    w_gene2.completions = genes
    if w_gene1.value not in genes:
        if 'APOE' in genes:
            w_gene1.value = 'APOE'
        else:
            w_gene1.value = genes[0]
    if w_gene2.value not in genes:
        if 'APOE' in genes:
            w_gene2.value = 'APOE'
        else:
            w_gene2.value = genes[0]


def update_facets():
    """Update interface for a specific dataset."""
    facets = get_facets()
    w_facet.options = facets
    if w_facet.value not in facets:
        # set
        w_facet.value = \
            [f for f in facets
             if not f.startswith('_')][0]


update_facets()
update_genes()


def get_data() -> pd.DataFrame:
    """Retrieve data from a dataset, gene & facet."""
    dataset_id = w_dataset_id.value
    gene1 = w_gene1.value
    gene2 = w_gene2.value
    facet = w_facet.value
    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene1}")

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene2}")

    data = pd.DataFrame(dict(
        gene1 = expset.get_gene(dataset_id, gene1)[:,0],
        gene2 = expset.get_gene(dataset_id, gene2)[:,0],
        obs = expset.get_meta(dataset_id, facet)[:,0],
        ))

    return data


def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]

def get_unique_obs(data):
    unique_obs = list(pd.DataFrame(data)['obs'].unique())
    return unique_obs



#
# Create plot
#

plot = figure()
source = ColumnDataSource(get_data())

if pd.DataFrame(source.data)['obs'].dtype == int:
    new_obs_column = expset.get_meta(w_dataset_id.value,w_facet.value,nobins=8)
    source.data['obs'] = new_obs_column[w_facet.value].to_numpy()
else:
    pass

index_cmap = factor_cmap('obs', Category10[len(list(pd.DataFrame(source.data)["obs"].unique()))], list(pd.DataFrame(source.data)["obs"].unique()))
glyph = plot.scatter(x='gene1', y='gene2', source=source,  legend_field="obs",
        fill_alpha=0.7, size=5, fill_color = index_cmap,width=0)

def cb_update_plot(attr, old, new):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, source, index_cmap,glyph
    data = get_data()
    dataset_id, dataset = get_dataset()
    facet = w_facet.value
    gene1 = w_gene1.value
    gene2 = w_gene2.value
    source.data = data

    if pd.DataFrame(source.data)['obs'].dtype == int:
        new_obs_column = expset.get_meta(w_dataset_id.value,w_facet.value,nobins=8)
        source.data['obs'] = new_obs_column[w_facet.value].to_numpy()
    else:
        pass

    unique_obs = get_unique_obs(data)
    index_cmap = factor_cmap('obs', Category10[len(unique_obs)], unique_obs)
    glyph.glyph.fill_color = index_cmap

   

    w_div_title_author.text = \
        f"""
        <ul>
          <li><b>Title:</b> {dataset['title']}</li>
          <li><b>Author:</b> {dataset['author']}</li>
          <li><b>Organism / Datatype:</b>
              {dataset['organism']} / {dataset['datatype']}</li>
        </ul>
        """
    w_download_filename.text = f"exp_{dataset_id}_{facet}_{gene1}_{gene2}.tsv"

    title = dataset['title'][:60]
    plot.title.text = (f"{gene1} vs {gene2} - "
                       f"({dataset_id}) {title}...")

    plot.xaxis.axis_label = f"{gene1}"
    plot.yaxis.axis_label = f"{gene2}"

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
    update_plot()


def cb_sibling_change(attr, old, new):
    lg.debug("Sibling change: " + new)
    w_dataset_id.value = new

##TODO: fix download button. only downloads old data
def cb_download():
    data = get_data()
    source.data = data
    print(source.data)
    return CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")

w_gene1.on_change("value", cb_update_plot)
w_gene2.on_change("value", cb_update_plot)
w_sibling.on_change("value", cb_sibling_change)
w_dataset_id.on_change("value", cb_dataset_change)
w_facet.on_change("value", cb_update_plot)
w_download.js_on_click(cb_download())


#
# Build the document
#
curdoc().add_root(
    column([
        row([w_gene1,w_gene2, w_facet, w_sibling, w_download],
            sizing_mode='stretch_width'),
        row([w_div_title_author],
            sizing_mode='stretch_width'),
        row([w_gene_not_found],
            sizing_mode='stretch_width'),
        row([plot],
            sizing_mode='stretch_width'),
        row([w_dataset_id, w_download_filename],),
    ], sizing_mode='stretch_width')
)

