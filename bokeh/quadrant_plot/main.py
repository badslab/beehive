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
from os.path import dirname, join


lg = logging.getLogger('QuadrantPlot')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Quadrant Plot'

VIEW_NAME = "quadrant_plot"

create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets(view_name=VIEW_NAME)

w_div_title_author = Div(text="",width = 300)


dataset_options = [(k, "{short_title}, {short_author}, {datatype}".format(**v))
                   for k, v in datasets.items()]
DATASET_NUMBER = 0

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[DATASET_NUMBER][0],
                             visible=True, width = 300)

w_genes = create_widget("genes",MultiChoice,default = ["TREM2","APOE"],options = [], title="Select genes to highlight",value_type = list,  height_policy = "fit")

w_category1 = create_widget("category1",Select, 
                    options=[], title="Select First Category", width = 150, height = 50)

w_category2 = create_widget("category2",Select,
                    options=[], title="Select Second Category", width = 150, height = 50)

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

    #filter and keep only lfcs..
    unique_vars = list(set([x.replace("__lfc","").replace("__lcpm","").replace("__cell_frac","") for x in vars]))
    print(unique_vars)
    print("****")
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
    # data["highlight"] = np.where(   (data["px"] < w_spinner.value) & (data["px"] != 0.0), "significant on x axis", data["highlight"])
    # data["highlight"] = np.where(   (data["py"] < w_spinner.value) & (data["py"] != 0.0), "significant on y axis", data["highlight"])
    # data["highlight"] = np.where((data["px"] < w_spinner.value) & (data["py"] < w_spinner.value)
    #                             & (data["px"] != 0.0) & (data["py"] != 0.0), 
    #                             "significant on both axes", data["highlight"])

    data["highlight"] = np.where(   (data["px"] < w_spinner.value), "significant on x axis", data["highlight"])
    data["highlight"] = np.where(   (data["py"] < w_spinner.value), "significant on y axis", data["highlight"])
    data["highlight"] = np.where((data["px"] < w_spinner.value) & (data["py"] < w_spinner.value), "significant on both axes", data["highlight"])

    data["highlight"] = data.apply(lambda x: highlight_genes(x),axis = 1)
    data["color"] = data.apply(lambda x: color_genes(x),axis = 1)
    data["size"] = data.apply(lambda x: size_genes(x),axis = 1)
    data.sort_values("color",inplace=True)
    return data


def color_genes(x):
    if x["highlight"] == "highlighted genes":
        return "red"
    elif x["highlight"] == "significant on x axis":
        return "orange"
    elif x["highlight"] == "significant on y axis":
        return "green"
    elif x["highlight"] == "significant on both axes":
        return "cyan"
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

def set_defaults():
    defaults_dict = expset.get_defaults(w_dataset_id.value,VIEW_NAME)
    if defaults_dict == {}:
        return

    for def_vals in defaults_dict[VIEW_NAME]:
        if def_vals.get("dataset"):
            default_dsid = def_vals.get("dataset")
            for i in range(len(w_dataset_id.options)):
                if default_dsid == w_dataset_id.options[i][0]:
                    w_dataset_id.value = w_dataset_id.options[i][0]
    
    update_vars()
    update_genes()
    for def_vals in defaults_dict[VIEW_NAME]:
        if def_vals.get("category1"):
            w_category1.value =  def_vals.get("category1")
        if def_vals.get("category2"):
            w_category2.value =  def_vals.get("category2")
        if def_vals.get("genes"):
            w_genes.value =  def_vals.get("genes")
        if def_vals.get("spinner"):
            w_spinner.value =  float(def_vals.get("spinner"))
    return

if curdoc().session_context.request.arguments == {}:
    set_defaults()

TOOLTIPS = [
            ('Log Fold Change on x-axis', '@x'),
            ('Log Fold Change on y-axis', '@y'),
            ('Adjusted p-value on x-axis', '@px'),
            ('Adjusted p-value on y-axis', '@py'),
            ('Gene:', '@gene'),
           ]

plot = figure(height = 400,output_backend = "webgl")
plot.add_tools(HoverTool(tooltips=TOOLTIPS))

x_diag = [x for x in range(-10,10)]
y_diag = x_diag

plot.line(x_diag,y_diag,line_dash = "dashed",color = "black")


data = get_data()
dataset_id, dataset = get_dataset()
categ1 = w_category1.value
categ2 = w_category2.value
w_download_filename.text = f"exp_{dataset_id}_{categ1}_{categ2}.csv"

#x-axis, y-axis, both, user, none

source = ColumnDataSource(data)

# COLOR_PALETTE = ["red","green","blue","purple", "grey"]
# plot.scatter(x = "x", y = "y", source = source, 
# color=factor_cmap('highlight', palette = COLOR_PALETTE, factors = HIGHLIGHTED_OPTIONS),
# legend_group = "highlight")

plot.scatter(x = "x", y = "y", source = source, 
color="color",
size = "size",
legend_field = "highlight")

plot.xaxis.axis_label = categ1
plot.yaxis.axis_label = categ2

v_func_text_size  = """
var new_xs = new Array(xs.length)
for(var i = 0; i < xs.length; i++) {
    new_xs[i] = text_sizes_map[xs[i]]
}
return new_xs
"""
HIGHLIGHTED_OPTIONS = ["significant on x axis","significant on y axis","significant on both axes","highlighted genes","other genes"]
TEXT_SIZES = ["0px","0px","0px","15px","0px"]

text_size_map = dict(zip(HIGHLIGHTED_OPTIONS,TEXT_SIZES))

categorical_text_size_transformer = CustomJSTransform(args={"text_sizes_map": text_size_map}, v_func = v_func_text_size)

labels = LabelSet(x='x', y='y', text='gene',
              x_offset=0, y_offset=5, source=source,
              text_font_size = transform('highlight',categorical_text_size_transformer),
              text_font_style = "bold")

plot.add_layout(labels)
def cb_update_plot(attr, old, new,type_change = None):
    curdoc().hold()
    global plot, source, data
    #changed p-value, no need to get new data
    if type_change == "p_sig":
        data["highlight"] = np.where(   (data["px"] < w_spinner.value), "significant on x axis", data["highlight"])
        data["highlight"] = np.where(   (data["py"] < w_spinner.value), "significant on y axis", data["highlight"])
        data["highlight"] = np.where((data["px"] < w_spinner.value) & (data["py"] < w_spinner.value), "significant on both axes", data["highlight"])

        data["highlight"] = data.apply(lambda x: highlight_genes(x),axis = 1)
        data["color"] = data.apply(lambda x: color_genes(x),axis = 1)
        data["size"] = data.apply(lambda x: size_genes(x),axis = 1)
        curdoc().unhold()
        return
    #rehighlighted some genes, no need to get data
    elif type_change == "rehighlight_genes":
        data["highlight"] = data.apply(lambda x: highlight_genes(x),axis = 1)
        data["color"] = data.apply(lambda x: color_genes(x),axis = 1)
        data["size"] = data.apply(lambda x: size_genes(x),axis = 1)
        curdoc().unhold()
        return
    #something else..
    else:
        data = get_data()
        source.data = data
        categ1 = w_category1.value
        categ2 = w_category2.value

        w_div_title_author.text = \
            f"""
            <ul>
            <li><b>Title:</b> {dataset['title']}</li>
            <li><b>Author:</b> {dataset['author']}</li>
            <li><b>Organism / Datatype:</b>
                {dataset['organism']} / {dataset['datatype']}</li>
            </ul>
            """

        plot.xaxis.axis_label = categ1
        plot.yaxis.axis_label = categ2
        w_download_filename.text = f"exp_{dataset_id}_{categ1}_{categ2}.csv"
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
w_genes.on_change("value",partial(cb_update_plot, type_change="rehighlight_genes"))
w_spinner.on_change("value",partial(cb_update_plot, type_change = "p_sig"))

w_dataset_id.on_change("value",cb_update_plot_new_dataset)


w_download.js_on_event("button_click", CustomJS(args=dict(source=source, file_name = w_download_filename),
                            code=open(join(dirname(__file__), "templates/download_general.js")).read()))


curdoc().add_root(row([
        column([
            row([w_genes, w_spinner]),
            row([w_category1,w_category2]),
            row([w_download,w_dataset_id]),
            row([w_div_title_author]),
                ]),
        column([plot],sizing_mode="stretch_width")
        ])
)