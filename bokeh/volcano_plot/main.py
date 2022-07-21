import logging
from functools import partial
import logging
from pprint import pprint
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource,MultiChoice
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, TextInput, Div,
                                  Button)
from bokeh.plotting import figure, curdoc


from beehive import config, util, expset

lg = logging.getLogger('ScatterExpression')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Volcano Plot'


create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets()

args = curdoc().session_context.request.arguments


w_div_title_author = Div(text="")


dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

DATASET_NUMBER = 9

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[DATASET_NUMBER][0],
                             visible=False,)

# Possible siblings of this dataset
siblings = expset.get_dataset_siblings(w_dataset_id.value)
sibling_options = []
for k, v in siblings.items():
    sname = f"{v['organism']} / {v['datatype']}"
    sibling_options.append((k, sname))

w_sibling = create_widget("view", Select,
                          options=sibling_options,
                          default=w_dataset_id.value,
                          update_url=False)

siblings = expset.get_dataset_siblings(w_dataset_id.value)


#make widgets for the following:

#select different category to view padj vs lfc
w_category = create_widget("category",Select, 
                    options=[], title="Select Category")

#text input to highlight select genes
w_text_input = TextInput(value= "", title="Gene names to be highlighted:")

#download button..
w_download = Button(label='Download', align='end')

w_download_filename = Div(text="", visible=False,
                          name="download_filename")

w_genes = MultiChoice(value=["TREM2", "APOE"], options=[],title="Select genes to highlight")


def get_genes():
    """Get available genes for a dataset."""
    dataset_id = w_dataset_id.value
    genes = sorted(list(expset.get_genes(dataset_id)))
    return genes

def update_genes():
    genes = get_genes()
    genes_options = [(x,x) for x in genes]
    w_genes.options = genes_options

def update_vars():
    """Update interface for a specific dataset."""
    vars = expset.get_varfields(w_dataset_id.value)
    unique_vars = list(set([x.replace("__lfc","").replace("__padj","") for x in vars]))
    #print(unique_vars)
    vars_options = [(x,x) for x in unique_vars]
    options = vars_options
    w_category.options = options
    if w_category.value not in [x[0] for x in options]:
        # set a default
        w_category.value = options[0][0]

update_vars()
update_genes()

def get_data() -> pd.DataFrame:
    dataset_id = w_dataset_id.value
    categ = w_category.value
    lg.warning(f"!! Getting data for {dataset_id} {categ}")

    data = expset.get_dedata_new(dataset_id,categ)
    data.rename(columns = {categ+"__lfc":"lfc", categ+"__padj":"padj"},inplace = True)

    return data

def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]


dataset_id, dataset = get_dataset()

plot = figure(height = 400)

data = get_data()
data["padj"] = np.log10(data["padj"])*-1
source = ColumnDataSource(data)

plot.scatter(x = "lfc", y = "padj", source = source)


def cb_update_plot(attr, old, new):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, source

    data = get_data()
    data["padj"] = np.log10(data["padj"])*-1
    source.data = data
    curdoc().unhold()


# convenience shortcut
update_plot = partial(cb_update_plot, attr=None, old=None, new=None)

# run it directly to ensure there are initial values
update_plot()

w_category.on_change("value", cb_update_plot)
w_genes.on_change("value",cb_update_plot)


def cb_dataset_change(attr, old, new):
    """Dataset change."""
    lg.info("Dataset Change to:" + new)
    curdoc().hold()
    update_vars()
    #update_genes()
    #update_plot()


def cb_sibling_change(attr, old, new):
    lg.debug("Sibling change: " + new)
    w_dataset_id.value = new


cb_download = CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")

w_download.js_on_click(cb_download)

curdoc().add_root(
    column([
        row([
            column([
                row([w_category,w_genes], 
                    sizing_mode='stretch_width'),
                row([w_sibling],
                    sizing_mode='stretch_width'),
                row([w_download],
                    sizing_mode='stretch_width')],   
                # row([w_download,w_x_axis,w_y_axis,w_regression],
                #     sizing_mode='stretch_width')],   
                sizing_mode='stretch_width')],
            sizing_mode='stretch_width'),
        row([w_div_title_author],
            sizing_mode='stretch_width'),
        row([plot],
            sizing_mode='stretch_width'),
        row([w_dataset_id, w_download_filename],),
    ], sizing_mode='stretch_width')
)





