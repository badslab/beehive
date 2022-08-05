import logging
from functools import partial
import logging
from venv import create
import pandas as pd
import numpy as np
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, RadioGroup, CheckboxGroup, Slider,MultiSelect, 
                            LinearColorMapper, ColorBar, Label,LinearAxis, Range1d)
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, Div,
                                  Button, AutocompleteInput)
from bokeh.plotting import figure, curdoc
from bokeh.models.transforms import CustomJSTransform
from bokeh.transform import transform

from beehive import config, util, expset
from bokeh.util.hex import cartesian_to_axial
from bokeh.palettes import Viridis256,Magma256

lg = logging.getLogger('ScatterExpression')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Hexbin Expression'

GENE_OPTION = 0
FACET_OPTION = 1
COLOR_POSSIBILITIES = ["None", "Density (counts)","Metadata (Gene or Numerical Facet)"]
LABELS_AXIS = ["Gene", "Numerical Facet"]
create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets()

args = curdoc().session_context.request.arguments


w_div_title_author = Div(text="")

#datasets with titles
dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

# Dataset
dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]


w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[0][0],
                             visible=False,)

#making widgets for the following:

#X and Y axes
w_gene1 = create_widget("geneX", AutocompleteInput,
                       completions=[], default='APOE', title = "Gene", case_sensitive=False)
w_gene2 = create_widget("geneY", AutocompleteInput,
                       completions=[], default='TREM2', title = "Gene", case_sensitive=False)
w_facet_numerical_1 = create_widget("num_facetX", Select,
                                    options=[], title="Numerical Facet")
w_facet_numerical_2 = create_widget("num_facetY", Select,
                                    options=[], title="Numerical Facet")

w_x_axis_radio = create_widget(
    "x_axis_radio", RadioGroup, labels=LABELS_AXIS, default=0, title="X-Axis", value_type=int)
w_y_axis_radio = create_widget(
    "y_axis_radio", RadioGroup, labels=LABELS_AXIS, default=0, title="Y-Axis", value_type=int)

widget_axes = [w_gene1,w_gene2, w_facet_numerical_1,w_facet_numerical_2]

#Z axis
LABELS_GROUPING = ["Gene", "Numerical Facet"]

#color by gene expression
w_gene3 = create_widget("geneZ", AutocompleteInput,
                       completions=[], default='TREM2', title = "Gene", case_sensitive=False)
#color by continuous facet
w_facet_numerical_3 = create_widget("num_facetZ", Select,
                                    options=[], title="Numerical Facet")
#color/ yes or no?
w_color_check = create_widget(
    "w_color_check", RadioGroup, labels=COLOR_POSSIBILITIES, default=0, title="Color By:", value_type=int)
#color by gene or numerical facet?
w_z_axis_radio = create_widget(
    "z_axis_radio", RadioGroup, labels=LABELS_AXIS, default=0, title="Grouped:", value_type=int)

#user selecting size
w_size_slider = create_widget("size_picker",Slider,start=0.01, end=1, default=0.1,step=0.01,title = "Hex Tile Size",value_type = float)

#subset selection of categorical obs
w_subset_select = create_widget("subset_categories",MultiSelect,default = [],value_type = list,options = [])

w_facet = create_widget("facet", Select, options=[],
                        title="Subset by Categories")


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

def update_genes():
    """Update genes widget for a dataset."""
    genes = get_genes()
    w_gene1.completions = genes
    w_gene2.completions = genes
    w_gene3.completions = genes
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
    if w_gene3.value not in genes:
        if 'APOE' in genes:
            w_gene3.value = 'APOE'
        else:
            w_gene3.value = genes[0]

update_genes()

def update_facets():
    """Update interface for a specific dataset."""
    options = expset.get_facet_options(w_dataset_id.value)
    w_facet.options = options
    if w_facet.value not in [x[0] for x in options]:
        # set a default
        w_facet.value = options[0][0]

update_facets()


def update_numerical_facets():
    """"Get and update numerical facets only. Helpful for the w_numerical_facet widget"""
    options = expset.get_facet_options_numerical(w_dataset_id.value)
    w_facet_numerical_1.options = options
    w_facet_numerical_2.options = options
    w_facet_numerical_3.options = options

    # some datasets might not have any numerical facets
    if not(options):
        return

    if w_facet_numerical_1.value not in [x[0] for x in options]:
        # set a default
        w_facet_numerical_1.value = options[0][0]
    if w_facet_numerical_2.value not in [x[0] for x in options]:
        # set a default
        w_facet_numerical_2.value = options[0][0]
    if w_facet_numerical_3.value not in [x[0] for x in options]:
        # set a default
        w_facet_numerical_3.value = options[0][0]

update_numerical_facets()

unique_obs = []
def get_data() -> pd.DataFrame:
    """Retrieve data from a dataset, gene & facet."""
    dataset_id = w_dataset_id.value
    gene1 = w_gene1.value
    gene2 = w_gene2.value
    gene3 = w_gene3.value
    num_facet1 = w_facet_numerical_1.value
    num_facet2 = w_facet_numerical_2.value
    num_facet3 = w_facet_numerical_3.value

    facet = w_facet.value

    geneX = None
    geneY = None
    geneZ = None
    num_facetX = None
    num_facetY = None
    num_facetZ = None


    x_axis = w_x_axis_radio.active
    y_axis = w_y_axis_radio.active
    z_axis = w_z_axis_radio.active

    if x_axis == GENE_OPTION:
        geneX = expset.get_gene(dataset_id, gene1)[:, 0]
    else:
        num_facetX = expset.get_meta(dataset_id, num_facet1, raw=True)[:, 0]
    if y_axis == GENE_OPTION:
        geneY = expset.get_gene(dataset_id, gene2)[:, 0]
    else:
        num_facetY = expset.get_meta(dataset_id, num_facet2, raw=True)[:, 0]

    geneZ = expset.get_gene(dataset_id, gene3)[:, 0]
    num_facetZ = expset.get_meta(dataset_id, num_facet3, raw=True)[:, 0]

    ##TODO maybe check for numerical/categorical here?
    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene1}")

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene2}")

    data = pd.DataFrame(dict(
        geneX = geneX,
        geneY = geneY,
        geneZ = geneZ,
        obs=expset.get_meta(dataset_id, facet)[:, 0],
        num_facetX = num_facetX,
        num_facetY = num_facetY,
        num_facetZ = num_facetZ,
        ))

    return data

def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]

def get_unique_obs(data):
    """Fetch the unique facets from the data"""
    global w_subset_select
    unique_obs = pd.DataFrame(data)['obs'].unique()
    w_subset_select.options = unique_obs.tolist()

def modify_data():
    data = get_data()
 
    global unique_obs, X_AXIS, Y_AXIS
    unique_obs = get_unique_obs(data)
    categories = w_subset_select.value
    SIZE = w_size_slider.value

    q_full, r_full = cartesian_to_axial(data[X_AXIS], data[Y_AXIS], size = SIZE, orientation = "pointytop")

    if len(categories) == 0:
        # q, r = cartesian_to_axial(data[X_AXIS], data[Y_AXIS], size = SIZE, orientation = "pointytop")
        q,r = q_full,r_full
    else:
        #filter out uneeded groups:
        data = data.loc[np.where(data["obs"].isin(categories))].reset_index(drop = True)
        q, r = cartesian_to_axial(data[X_AXIS], data[Y_AXIS], size = SIZE, orientation = "pointytop")
    df = pd.DataFrame(dict(r=r,q=q,avg_exp = None))

    groups = df.groupby(["q","r"])
    dicts = []
    ###assign axial coordinates to get average gene expression###
    #need to che

    for key,val in groups.groups.items():
        q,r = key
        counts = len(val)

        coloring_scheme = data[Z_AXIS].loc[groups.groups.get((q,r))].mean()
        dicts = dicts + [{"q" : q,"r" : r,"counts" : counts, "coloring_scheme" :coloring_scheme,"index_list":val}]

    final_result = pd.DataFrame(dicts)
    data_full = pd.DataFrame(dict(r=r_full,q=q_full))
    return final_result, data_full


def get_mapper(feature,data,palette):
    return  LinearColorMapper(
    palette=palette,
    low=np.percentile(np.array(data[feature]),1),
    high=np.percentile(np.array(data[feature]),99))

#
# Create plot
#
X_AXIS = "geneX" if w_x_axis_radio.active == GENE_OPTION else "num_facetX"
Y_AXIS = "geneY" if w_y_axis_radio.active == GENE_OPTION else "num_facetY"
Z_AXIS = "geneZ" if w_z_axis_radio.active == GENE_OPTION else "num_facetZ"


plot = figure(output_backend = "webgl")
dataset_id, dataset = get_dataset()
data,hex_full = modify_data()

mapper_metadata = get_mapper("coloring_scheme", data, [x for x in Magma256][::-1])
mapper_counts = get_mapper("counts",data,[x for x in Viridis256][::-1])

source = ColumnDataSource(data)
source_full = ColumnDataSource(hex_full)
ALPHA = 0.8
SIZE = w_size_slider.value

length = len(sorted(np.unique(source.data["counts"]).tolist()))

plot.hex_tile(q="q",r="r",size = SIZE,source=source_full,
    color  = "#D3D3D3",line_color = None, alpha = ALPHA) #grey plot all

visibility_counts = True
visibility_metadata = True

if COLOR_POSSIBILITIES[w_color_check.active] == "None":
    value_text = ""
    tick_min_text = ""
    tick_max_text = ""
    visibility_counts = False
    visibility_metadata = False
    pass
    
elif COLOR_POSSIBILITIES[w_color_check.active] == "Density (counts)":
    plot.hex_tile(q="q",r="r",size = SIZE,source=source,
    line_color = None, 
    color   = {'field': 'counts', 'transform': mapper_counts} , alpha = ALPHA)
    value_text = "counts"
    tick_min_text = f'{np.round(np.percentile(np.array(data["counts"]),1),2)}'
    tick_max_text = f'{np.round(np.percentile(np.array(data["counts"]),99),2)}'
    visibility_counts = True
    visibility_metadata = False


else:
    plot.hex_tile(q="q",r="r",size = SIZE,source=source,
    line_color = None, 
    color   = {'field': 'coloring_scheme', 'transform': mapper_metadata}, alpha = ALPHA)
    value_text = w_gene3.value if w_z_axis_radio.active == GENE_OPTION else w_facet_numerical_3.value
    tick_min_text = f'{np.round(np.percentile(np.array(data["coloring_scheme"]),1),2)}'
    tick_max_text = f'{np.round(np.percentile(np.array(data["coloring_scheme"]),99),2)}'
    visibility_counts = False
    visibility_metadata = True


#1. colorbar
color_bar_meta = ColorBar(color_mapper=mapper_metadata, location=(
    0, 0), major_label_text_font_size="0px", major_tick_in=2)
color_bar_counts = ColorBar(color_mapper=mapper_counts, location=(
    0, 0), major_label_text_font_size="0px", major_tick_in=2)
# 2. colorbar title

colorbar_text = Label(text=value_text, x=-20, y=520, x_units='screen',
                      y_units='screen', text_align="center", text_font_style="italic", text_font_size="12px")
# 3. min value
tick_min = Label(text=tick_min_text, x=-40, y=15, y_units='screen',
                 x_units="screen", text_align="left", text_baseline="middle", text_font_size="12px")
# 4. max value
tick_max = Label(text=tick_max_text, x=-50, y=510, y_units='screen',
                 x_units="screen", text_align="left", text_baseline="middle", text_font_size="12px")

plot.add_layout(color_bar_meta, "right")
plot.add_layout(color_bar_counts, "right")
plot.add_layout(colorbar_text, "right")
plot.add_layout(tick_min, "right")
plot.add_layout(tick_max, "right")

color_bar_meta.visible = visibility_metadata
color_bar_counts.visible = visibility_metadata

if COLOR_POSSIBILITIES[w_color_check.active] == "None":
    colorbar_text.visible = False
    tick_min.visible = False
    tick_max.visible = False

x_label = ""
y_label = ""

def cb_update_plot(attr, old, new,type_change):
    """Populate and update the plot."""
    curdoc().hold()
    global plot,source, X_AXIS, Y_AXIS, Z_AXIS, SIZE, ALPHA, widget_axes,x_label,y_label,color_bar, colorbar_text, tick_max,tick_min


    X_AXIS = "geneX" if w_x_axis_radio.active == GENE_OPTION else "num_facetX"
    Y_AXIS = "geneY" if w_y_axis_radio.active == GENE_OPTION else "num_facetY"
    Z_AXIS = "geneZ" if w_z_axis_radio.active == GENE_OPTION else "num_facetZ"


    data,hex_full = modify_data()
    # get_unique_obs(data)
    dataset_id, dataset = get_dataset()

    facet = w_facet.value
    source = ColumnDataSource(data)
    source_full = ColumnDataSource(hex_full)
    SIZE = w_size_slider.value

    color_bar_counts.visible = False
    color_bar_meta.visible = False
    colorbar_text.visible = False
    tick_min.visible = False
    tick_max.visible = False

    #update plot now
    plot.renderers = []

    mapper_metadata = get_mapper("coloring_scheme", data, [x for x in Magma256][::-1])
    mapper_counts = get_mapper("counts",data,[x for x in Viridis256][::-1])

    plot.hex_tile(q="q",r="r",size = SIZE,source=source_full,
        color  = "#D3D3D3",line_color = None, alpha = ALPHA) #grey plot all


    if COLOR_POSSIBILITIES[w_color_check.active] == "None":
        pass
    else:
        if COLOR_POSSIBILITIES[w_color_check.active] == "Density (counts)":
            plot.hex_tile(q="q",r="r",size = SIZE,source=source,
            line_color = None, 
            color   = {'field': 'counts', 'transform': mapper_counts}, alpha = ALPHA)
            colorbar_text.text = "counts"
            tick_min.text = f'{np.round(np.percentile(np.array(data["counts"]),1),2)}'
            tick_max.text = f'{np.round(np.percentile(np.array(data["counts"]),99),2)}'
            color_bar_counts.visible = True

        else:
            plot.hex_tile(q="q",r="r",size = SIZE,source=source,
            line_color = None, 
            color   = {'field': 'coloring_scheme', 'transform': mapper_metadata}, alpha = ALPHA)
            colorbar_text.text = w_gene3.value if w_z_axis_radio.active == GENE_OPTION else w_facet_numerical_3.value
            tick_min.text = f'{np.round(np.percentile(np.array(data["coloring_scheme"]),1),2)}'
            tick_max.text = f'{np.round(np.percentile(np.array(data["coloring_scheme"]),99),2)}'
            color_bar_meta.visible = True

        # make them visible again
        colorbar_text.visible = True
        tick_min.visible = True
        tick_max.visible = True

    w_div_title_author.text = \
        f"""
        <ul>
          <li><b>Title:</b> {dataset['title']}</li>
          <li><b>Author:</b> {dataset['author']}</li>
          <li><b>Organism / Datatype:</b>
              {dataset['organism']} / {dataset['datatype']}</li>
        </ul>
        """

    title = dataset['title'][:60]

    for ind,widget in enumerate(widget_axes):
        if X_AXIS in widget.name:
            x_label = widget_axes[ind].value
        if Y_AXIS in widget.name:
            y_label = widget_axes[ind].value

    plot.title.text = (f"{x_label} vs {y_label} - "
                       f"({dataset_id}) {title}...")

    #only change main x and y axes labels...
    plot.xaxis[0].axis_label = f"{x_label}"
    plot.yaxis[0].axis_label = f"{y_label}"

    w_download_filename.text = f"exp_{dataset_id}_{facet}_{x_label}_{y_label}.tsv"
    
    curdoc().unhold()


# convenience shortcut
update_plot = partial(cb_update_plot, attr=None, old=None, new=None,type_change=None)

# run it directly to ensure there are initial values
update_plot()


def cb_dataset_change(attr, old, new):
    """Dataset change."""
    lg.info("Dataset Change to:" + new)
    curdoc().hold()
    update_facets()
    update_genes()
    update_plot()

cb_download = CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")

w_gene1.on_change("value",partial(cb_update_plot,type_change = None))
w_gene2.on_change("value", partial(cb_update_plot,type_change=None))
w_gene3.on_change("value",partial(cb_update_plot,type_change=None))
w_facet.on_change("value", partial(cb_update_plot,type_change=None))
w_facet_numerical_1.on_change("value", partial(cb_update_plot,type_change=None))
w_facet_numerical_2.on_change("value", partial(cb_update_plot,type_change=None))
w_facet_numerical_3.on_change("value", partial(cb_update_plot,type_change=None))
w_subset_select.on_change("value",partial(cb_update_plot,type_change = None))
w_size_slider.on_change("value",partial(cb_update_plot,type_change = None))
w_x_axis_radio.on_change("active",partial(cb_update_plot,type_change = None))
w_y_axis_radio.on_change("active",partial(cb_update_plot,type_change = None))
w_z_axis_radio.on_change("active",partial(cb_update_plot,type_change = None))
w_color_check.on_change("active",partial(cb_update_plot,type_change = None))
w_dataset_id.on_change("value", cb_dataset_change)
w_download.js_on_click(cb_download)

curdoc().add_root(
    column([
        row([
            column([
                row([w_gene1,w_gene2,w_gene3], 
                    sizing_mode='stretch_width'),
                row([w_facet_numerical_1,w_facet_numerical_2,w_facet_numerical_3],
                    sizing_mode='stretch_width'),
                row([w_x_axis_radio,w_y_axis_radio,w_z_axis_radio],
                    sizing_mode='stretch_width'),
                row([w_facet,w_subset_select,w_size_slider, w_download],
                    sizing_mode='stretch_both'),
                row([w_color_check],
                    sizing_mode='stretch_width')],   
                sizing_mode='stretch_width')],
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