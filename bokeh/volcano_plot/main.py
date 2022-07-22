import logging
from functools import partial
import logging
from pprint import pprint
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource,MultiChoice, HoverTool, Spinner, Range1d
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, TextInput, Div,
                                  Button)
from bokeh.plotting import figure, curdoc
from bokeh.transform import factor_cmap, transform
from bokeh.models.transforms import CustomJSTransform
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


#download button..
w_download = Button(label='Download', align='end')

w_download_filename = Div(text="", visible=False,
                          name="download_filename")


#w_genes = MultiChoice(value=["TREM2", "APOE"], options=[],title="Select genes to highlight")
w_genes = create_widget("genes",MultiChoice,default = ["TREM2","APOE"],options = [], title="Select genes to highlight",value_type = list)
w_x_range = create_widget("x_range",Spinner,title = "Select Min & Max Absolute X-range",low=0, high = 4,step=0.1,default=2,value_type = float)
w_y_range = create_widget("y_range",Spinner,title = "Select Max Y-range",low=0, high = 300,step=1,default=100,value_type = float)


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
    #drop last column,, the one for genes.
    vars = vars[:-1]
    unique_vars = list(set([x.replace("__lfc","").replace("__padj","") for x in vars]))
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
    data.columns = ["lfc","padj","gene"]
    return data


def highlight_genes(x):
    if x["gene"] in w_genes.value:
        return "Yes"
    else:
        return "No"


def modify_data():
    data = get_data()

    data["padj"] = np.log10(data["padj"])*-1
    data["highlight"] = data.apply(lambda x: highlight_genes(x),axis = 1)

    #take care of x y ranges
    x_range = w_x_range.value
    y_range = w_y_range.value

    return data


def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]


dataset_id, dataset = get_dataset()
tooltips = [
            ('Log Fold Change:', '@lfc'),
            ('P-value adjsuted:', '@padj'),
            ('Gene:', '@gene'),
           ]



##Plot###
plot = figure(height = 400)
plot.add_tools(HoverTool(tooltips=tooltips))

data = modify_data()

source = ColumnDataSource(data)

v_func_alpha  = """
var new_xs = new Array(xs.length)
for(var i = 0; i < xs.length; i++) {
    new_xs[i] = alpha_map[xs[i]]
}
return new_xs
"""

v_func_size  = """
var new_xs = new Array(xs.length)
for(var i = 0; i < xs.length; i++) {
    new_xs[i] = size_map[xs[i]]
}
return new_xs
"""


alphas = [1,0.01]
alpha_map = dict(zip(["Yes","No"], alphas))
categorical_alpha_transformer = CustomJSTransform(args={"alpha_map": alpha_map}, v_func=v_func_alpha)

sizes = [10,2]
size_map = dict(zip(["Yes","No"],sizes))
categorical_size_transformer =  CustomJSTransform(args={"size_map": size_map},v_func=v_func_size)

plot.scatter(x = "lfc", y = "padj", source = source,
color=factor_cmap('highlight', palette = ["red","blue"], factors = ["Yes","No"]),
fill_alpha = transform('highlight',categorical_alpha_transformer),
size = transform('highlight',categorical_size_transformer))

plot.update(x_range = Range1d(w_x_range.value*-1, w_x_range.value), y_range = Range1d(0, w_y_range.value))
def cb_update_plot(attr, old, new):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, source

    data = modify_data()

    source.data = data

    plot.x_range.update(start = w_x_range.value*-1, end = w_x_range.value)
    plot.y_range.update(start = 0, end = w_y_range.value)

    curdoc().unhold()


# convenience shortcut
update_plot = partial(cb_update_plot, attr=None, old=None, new=None)

# run it directly to ensure there are initial values
update_plot()

w_category.on_change("value", cb_update_plot)
w_genes.on_change("value",cb_update_plot)
w_x_range.on_change("value",cb_update_plot)
w_y_range.on_change("value",cb_update_plot)

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
                row([w_genes], 
                    sizing_mode='stretch_width'),
                row([w_category,w_sibling],
                    sizing_mode='stretch_width'),
                row([w_x_range,w_y_range,w_download],
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





