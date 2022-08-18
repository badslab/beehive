import logging
from functools import partial
import logging
from pprint import pprint
import numpy as np
import pandas as pd
import polars as pl

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource,MultiChoice, HoverTool, Spinner, Range1d, LabelSet
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, Div,
                                  Button)
from bokeh.plotting import figure, curdoc
from bokeh.transform import factor_cmap, transform
from bokeh.models.transforms import CustomJSTransform
from beehive import config, util, expset


lg = logging.getLogger('QuadrantPlot')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Quadrant Plot'

VIEW_NAME = "quadrant_plot"

create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets(view_name=VIEW_NAME)

args = curdoc().session_context.request.arguments

w_div_title_author = Div(text="",width = 100)

PVAL_CUTOFF = 0.0001

dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

DATASET_NUMBER = 0

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[DATASET_NUMBER][0],
                             visible=True, width = 300)

w_genes = create_widget("genes",MultiChoice,default = ["TREM2","APOE"],options = [], title="Select genes to highlight",value_type = list,  
height = 100, width = 300)

w_category1 = create_widget("category1",Select, 
                    options=[], title="Select First Category", width = 150, height = 50)

w_category2 = create_widget("category2",Select,
                    options=[], title="Select Second Category", width = 150, height = 50)

w_download = Button(label='Download', align='end', width = 150, height = 50)
w_download_filename = Div(text="", visible=False,
                          name="download_filename")


def get_genes():
    """Get available genes for a dataset."""
    dataset_id = w_dataset_id.value
    genes = sorted(list(expset.get_genes(dataset_id)))
    return genes

def update_genes():
    """"fetch and update the genes of the widget"""
    genes = get_genes()
    genes_options = [(x,x) for x in genes]
    w_genes.options = genes_options

def update_vars():
    """Update interface for a specific dataset."""
    vars = expset.get_varfields(w_dataset_id.value)
    #drop last column,, the one for genes.
    vars = vars[:-1]
    vars = list(filter(lambda x: "__padj" not in x,vars))

    #filter and keep only lfcs..
    unique_vars = list(set([x.replace("__lfc","") for x in vars]))

    vars_options = [(x,x) for x in unique_vars]
    
    w_category1.options = vars_options
    w_category2.options = vars_options

    if w_category1.value not in [x[0] for x in vars_options]:
        w_category1.value = vars_options[0][0]

    if w_category2.value not in [x[0] for x in vars_options]:
        # set a default
        w_category2.value = vars_options[0][0]


update_vars()
update_genes()

def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]



def get_data() -> pd.DataFrame:
    """Get data of lfc on x axis and lfc on y axis and  gene names associated with that var"""
    dataset_id = w_dataset_id.value
    categ1 = w_category1.value
    categ2 = w_category2.value
    lg.warning(f"!! Getting data for {dataset_id} {categ1} {categ2}")

    data = expset.get_dedata_quadrant(dataset_id,categ1,categ2)
    data.columns = ["x","y","px","py","gene"]
    #dropping NaN values for now.
    data.dropna(inplace = True)


    data["highlight"] = np.where(True, "other genes","other genes")

    #exclude padj == 0?
    data["highlight"] = np.where(   (data["px"] < PVAL_CUTOFF) & (data["px"] != 0.0), "significant on x axis", data["highlight"])
    data["highlight"] = np.where(   (data["py"] < PVAL_CUTOFF) & (data["py"] != 0.0), "significant on y axis", data["highlight"])
    data["highlight"] = np.where((data["px"] < PVAL_CUTOFF) & (data["py"] < PVAL_CUTOFF)
                                & (data["px"] != 0.0) & (data["py"] != 0.0), 
                                "significant on both axes", data["highlight"])
    data["highlight"] = data.apply(lambda x: highlight_genes(x),axis = 1)
    data["color"] = data.apply(lambda x: color_genes(x),axis = 1)
    data["size"] = data.apply(lambda x: size_genes(x),axis = 1)
    return data


def color_genes(x):
    if x["highlight"] == "highlighted genes":
        return "purple"
    elif x["highlight"] == "significant on x axis":
        return "red"
    elif x["highlight"] == "significant on y axis":
        return "green"
    elif x["highlight"] == "significant on both axes":
        return "blue"
    else:
        return "gray"
       
def highlight_genes(x):
    """"Helper function to label genes as Yes or No in the highlight column"""
    if x["gene"] in w_genes.value:
        return "highlighted genes"
    else:
        return x["highlight"]

def size_genes(x):
    if x["highlight"] == "other genes":
        return 3
    else:
        return 6

TOOLTIPS = [
            ('Log Fold Change on x-axis', '@x'),
            ('Log Fold Change on y-axis', '@y'),
            ('Adjusted p-value on x-axis', '@px'),
            ('Adjusted p-value on y-axis', '@py'),
            ('Gene:', '@gene'),
           ]

plot = figure(height = 400,output_backend = "webgl")
plot.add_tools(HoverTool(tooltips=TOOLTIPS))

data = get_data()
dataset_id, dataset = get_dataset()
categ1 = w_category1.value
categ2 = w_category2.value

#x-axis, y-axis, both, user, none

source = ColumnDataSource(data)

# COLOR_PALETTE = ["red","green","blue","purple", "grey"]
# HIGHLIGHTED_OPTIONS = ["significant on x axis","significant on y axis","significant on both axes","highlighted genes","other genes"]
# plot.scatter(x = "x", y = "y", source = source, 
# color=factor_cmap('highlight', palette = COLOR_PALETTE, factors = HIGHLIGHTED_OPTIONS),
# legend_group = "highlight")

plot.scatter(x = "x", y = "y", source = source, 
color="color",
size = "size",
legend_field = "highlight")

def cb_update_plot(attr, old, new,type_change = None):
    curdoc().hold()
    global plot, source

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


    curdoc().unhold()


update_plot = partial(cb_update_plot, attr=None, old=None, new=None,type_change = None)

# run it directly to ensure there are initial values
update_plot()


def cb_update_plot_new_dataset(attr,old,new):
    curdoc().hold()
    
    global plot, source
    
    #update the w_category widget
    update_vars()
    #update the w_gene widget
    update_genes()
    #fetch data again
    data = get_data()

    source.data = data
    curdoc().unhold()


w_category1.on_change("value", partial(cb_update_plot))
w_category2.on_change("value", partial(cb_update_plot))

w_genes.on_change("value",partial(cb_update_plot))

w_dataset_id.on_change("value",cb_update_plot_new_dataset)

cb_download = CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")

w_download.js_on_click(cb_download)


curdoc().add_root(row([
        column([
            row([w_genes]),
            row([w_category1,w_category2]),
            row([w_download,w_dataset_id]),
            row([w_div_title_author],sizing_mode="scale_both"),
                ]),
        column([plot],sizing_mode="scale_both")
        ],sizing_mode="stretch_both")
)








