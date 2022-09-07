import logging
from functools import partial
import logging
import numpy as np
import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource,MultiChoice, HoverTool, Spinner, Range1d, LabelSet
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, Div,
                                  Button)
from bokeh.plotting import figure, curdoc
from beehive import config, util, expset
from bokeh.models.transforms import CustomJSTransform
from bokeh.transform import transform


lg = logging.getLogger('MeanAbundance')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Mean Abundance'


VIEW_NAME = "mean_abundance"

create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets(view_name=VIEW_NAME)

args = curdoc().session_context.request.arguments

w_div_title_author = Div(text="",width = 300)

dataset_options = [(k, "{short_title}, {short_author}, {datatype}".format(**v))
                   for k, v in datasets.items()]
DATASET_NUMBER = 0

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[DATASET_NUMBER][0],
                             visible=True, width = 300)

w_genes = create_widget("genes",MultiChoice,default = ["TREM2","APOE"],options = [], title="Select genes to highlight",value_type = list,  
height_policy = "fit")


w_category = create_widget("category",Select, 
                    options=[], title="Select Category", width = 150, height = 50)

w_download = Button(label='Download', align='end', width = 150, height = 50)
w_download_filename = Div(text="", visible=False,
                          name="download_filename")

w_spinner = create_widget("spinner",Spinner, title="p value significance", low = 0, high = 0.10, step = 0.0001, default = 0.05,value_type = float, width = 100)

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
    #keep only lfcs that have lcpm (abundance count)
    results = [string for string in vars if string.endswith("__lcpm")]
    unique_vars = list(set([x.replace("__lcpm","") for x in results]))
    vars_options = [(x,x) for x in unique_vars]
    
    w_category.options = vars_options

    if w_category.value not in [x[0] for x in vars_options]:
        w_category.value = vars_options[0][0]

update_vars()
update_genes()

def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]



def get_data() -> pd.DataFrame:
    """Get data of lfc on x axis and lfc on y axis and  gene names associated with that var"""
    dataset_id = w_dataset_id.value
    categ = w_category.value
    lg.warning(f"!! Getting data for {dataset_id} {categ}")

    data = expset.get_dedata_abundance(dataset_id,categ)

    data.columns = ["lfc","padj","abundance","gene"]
    #dropping NaN values for now.
    data.dropna(inplace = True)
    return data


TOOLTIPS = [
            ('Log Fold Change', '@lfc'),
            ('Adjusted p-value', '@padj'),
            ('Abundance', '@abundance'),
            ('Gene:', '@gene'),
           ]


plot = figure(output_backend = "webgl")
plot.add_tools(HoverTool(tooltips=TOOLTIPS))

data = get_data()
dataset_id, dataset = get_dataset()
categ = w_category.value


source = ColumnDataSource(data)
plot.scatter(x = "abundance", y = "lfc", source = source)

def cb_update_plot(attr, old, new,type_change = None):
    curdoc().hold()
    global plot, source, data
    #changed p-value, no need to get new data
    #rehighlighted some genes, no need to get data
    #something else..
    data = get_data()
    source.data = data
    categ = w_category.value

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


w_category.on_change("value", partial(cb_update_plot))
w_genes.on_change("value",partial(cb_update_plot))
w_spinner.on_change("value",partial(cb_update_plot))

w_dataset_id.on_change("value",cb_update_plot_new_dataset)

cb_download = CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")

w_download.js_on_click(cb_download)


curdoc().add_root(row([
        column([
            row([w_genes, w_spinner]),
            row([w_category,w_download]),
            row([w_dataset_id]),
            row([w_div_title_author]),
                ]),
        column([plot],sizing_mode="stretch_width")
        ])
)