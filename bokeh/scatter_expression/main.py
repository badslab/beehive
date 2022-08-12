from ctypes import sizeof
from enum import unique
import logging
from functools import partial
import logging
import pandas as pd
from scipy import stats
import numpy as np
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CheckboxGroup, RadioGroup, LinearColorMapper, ColorBar, Label,Slider)
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, Div,
                                  Button, AutocompleteInput)
from bokeh.plotting import figure, curdoc
from bokeh.transform import factor_cmap
from bokeh.palettes import Category20


from beehive import config, util, expset
# import traceback


lg = logging.getLogger('ScatterExpression')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Scatter Expression'

GENE_OPTION = 0
FACET_OPTION = 1

create_widget = partial(util.create_widget, curdoc=curdoc())

datasets = expset.get_datasets()

args = curdoc().session_context.request.arguments


# TODO
# clear cookies for user/password? ==> in auth not here..
# print(curdoc().session_context.request.cookies)

w_div_title_author = Div(text="")

# datasets with titles

# Dataset
dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

# TODO setting manually
DATASET_NUMBER = 0

# Widgets:
w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[DATASET_NUMBER][0],
                             visible=True,)

w_gene1 = create_widget("geneX", AutocompleteInput,
                        completions=[], default='APOE', title="Gene", case_sensitive=False,width=100)
w_gene2 = create_widget("geneY", AutocompleteInput,
                        completions=[], default='TREM2', title="Gene", case_sensitive=False,width=100)
w_facet_numerical_1 = create_widget("num_facetX", Select,
                                    options=[], title="Numerical Facet",width=100)
w_facet_numerical_2 = create_widget("num_facetY", Select,
                                    options=[], title="Numerical Facet",width=100)

widget_axes = [w_gene1, w_gene2, w_facet_numerical_1, w_facet_numerical_2]

w_regression = CheckboxGroup(labels=["Regression Lines"], active=[])

# Widget fixed options:
FIXED_OPTIONS = [("geneX", "Gene X"), ("geneY", "Gene Y"), ("num_facetX",
                                                            "Numerical Facet on X"), ("num_facetY", "Numerical Facet on Y")]

LABELS_AXIS = ["Gene", "Numerical Facet"]
LABELS_GROUPING = ["categorical facet", "gene","numerical facet"]

w_x_axis_radio = create_widget(
    "x_axis_radio", RadioGroup, labels=LABELS_AXIS, default=0, title="X-Axis", value_type=int,width=120)
w_y_axis_radio = create_widget(
    "y_axis_radio", RadioGroup, labels=LABELS_AXIS, default=0, title="Y-Axis", value_type=int,width=120)

# categorical facet grouping of data... or gene3 grouping
w_facet = create_widget("facet", Select, options=[],
                        title="Group by Categories",width=100)

w_gene3 = create_widget("geneZ", AutocompleteInput,
                        completions=[], default="APOE", title="Group by Gene Expression", case_sensitive=False,width=150)
                        
w_facet_numerical_3 = create_widget("num_facetZ", Select,
                                    options=[], title="Numerical Facet",width=100)

w_category_radio = create_widget(
    "category_radio", RadioGroup, labels=LABELS_GROUPING, default=0, title="Grouped:", value_type=int,width=100)

# download button..
w_download = Button(label='Download', align='end',width=100)

w_download_filename = Div(text="", visible=False,
                          name="download_filename")

w_alpha_slider = create_widget("alpha_picker",Slider,start=0.1, end=1, default=0.7,step=0.05,title = "Opacity",value_type = float,width=100)
w_size_slider = create_widget("size_picker",Slider,start=2, end=15, default=5,step=1,title = "Points Size",value_type = int,width=100)

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


def update_facets():
    """Update interface for a specific dataset."""
    options = expset.get_facet_options(w_dataset_id.value)
    w_facet.options = options
    if w_facet.value not in [x[0] for x in options]:
        # set a default
        w_facet.value = options[0][0]


update_facets()
update_genes()


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
    if len(options) == 0:
        w_category_radio.update(labels = ["categorical facet", "gene"])

update_numerical_facets()


def get_data() -> pd.DataFrame:
    """Retrieve data from a dataset, gene & facet."""
    dataset_id = w_dataset_id.value
    gene1 = w_gene1.value
    gene2 = w_gene2.value
    facet = w_facet.value
    num_facet1 = w_facet_numerical_1.value
    num_facet2 = w_facet_numerical_2.value
    num_facet3 = w_facet_numerical_3.value

    # all data will have 6 columns
    # initially set at None.
    geneX = None
    geneY = None
    num_facetX = None
    num_facetY = None
    geneZ = None
    gene3 = w_gene3.value
    num_facetZ = None

    # only get the data that we want to visualize
    # x_axis and y_axis values indicate which
    x_axis = w_x_axis_radio.active
    y_axis = w_y_axis_radio.active

    if x_axis == GENE_OPTION:
        geneX = expset.get_gene(dataset_id, gene1)[:, 0]
    else:
        num_facetX = expset.get_meta(dataset_id, num_facet1, raw=True)[:, 0]

    if y_axis == GENE_OPTION:
        geneY = expset.get_gene(dataset_id, gene2)[:, 0]
    else:
        num_facetY = expset.get_meta(dataset_id, num_facet2, raw=True)[:, 0]

    # always fetched => for the colorbar when coloring
    # with a gene expression
    if gene3:
        geneZ = expset.get_gene(dataset_id, gene3)[:, 0]
    
    if num_facet3:
        num_facetZ = expset.get_meta(dataset_id, num_facet3, raw=True)[:, 0]

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene1}")

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene2}")

    data = pd.DataFrame(dict(
        geneX=geneX,
        geneY=geneY,
        obs=expset.get_meta(dataset_id, facet)[:, 0],
        num_facetX=num_facetX,
        num_facetY=num_facetY,
        geneZ=geneZ,
        num_facetZ = num_facetZ
    ))
    if len(w_facet_numerical_1.options) == 0:
        data.drop("num_facetZ",axis=1,inplace = True)

    return data


def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]


def get_unique_obs(data):
    """Fetch the unique facets from the data"""
    unique_obs = pd.DataFrame(data)['obs'].unique()
    return unique_obs


def get_mapper():
    """"Set up palette for ColorBar"""

    if LABELS_GROUPING[w_category_radio.active] == "gene" or len(w_facet_numerical_1.options) == 0 :
        mapper = LinearColorMapper(
            palette='Magma256',
            low=data["geneZ"].min(),
            high=data["geneZ"].max())
    else:
        mapper = LinearColorMapper(
            palette='Magma256',
            low=data["num_facetZ"].min(),
            high=data["num_facetZ"].max())

    return mapper

#
# Create plot
#


plot = figure(output_backend="webgl",width=1000)

data = get_data()

unique_obs = get_unique_obs(data)
dataset_id, dataset = get_dataset()

# palette for the categorical obs
if len(unique_obs) < 3:
    palette = ["blue","red"]
else:
    palette = Category20[len(unique_obs)]

index_cmap = factor_cmap('obs', palette, unique_obs)

# initial set up of which x-axis and y-axis we plot
X_AXIS = "geneX"
Y_AXIS = "geneY"
Z_AXIS = LABELS_GROUPING[w_category_radio.active]

# for colorbar
mapper = get_mapper()
ALPHA = w_alpha_slider.value
SIZE = w_size_slider.value

# Are we plotting the coloring with obs categorical or with obs numerical?
if Z_AXIS == "categorical facet":
    # Plot multiple glyphs, each glyph is for 1 unique obs group
    for index, obs in enumerate(unique_obs):

        new_source = ColumnDataSource(data.loc[(data.obs == obs)])

        plot.scatter(x=X_AXIS, y=Y_AXIS, source=new_source,  legend_label=obs,
                     fill_alpha=ALPHA, size=SIZE, width=0, fill_color=index_cmap["transform"].palette[index])
    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"
elif Z_AXIS == "gene":
    plot.scatter(x=X_AXIS, y=Y_AXIS, source=ColumnDataSource(data),
                 fill_alpha=ALPHA, size=SIZE, width=0, fill_color={'field': 'geneZ', 'transform': mapper})
else:
    plot.scatter(x=X_AXIS, y=Y_AXIS, source=ColumnDataSource(data),
                 fill_alpha=ALPHA, size=ALPHA, width=0, fill_color={'field': 'num_facetZ', 'transform': mapper})

# Colorbar###:
# 4 components:
if Z_AXIS == "gene":
    text_label = w_gene3.value
    min_text = data["geneZ"].min()
    max_text = data["geneZ"].max()
elif Z_AXIS == "numerical facet":
    text_label = w_facet_numerical_3.value
    min_text = data["num_facetZ"].min()
    max_text = data["num_facetZ"].max()
else:
    text_label = ""
    min_text = "1"
    max_text = "1"
# 1.colorbar
color_bar = ColorBar(color_mapper=mapper, location=(
    0, 0), major_label_text_font_size="0px", major_tick_in=2)
plot.add_layout(color_bar, "right")
# 2. colorbar title
colorbar_text = Label(text=f'{text_label}', x=-20, y=520, x_units='screen',
                      y_units='screen', text_align="center", text_font_style="italic", text_font_size="12px")
plot.add_layout(colorbar_text, "right")
# 3. min value

tick_min = Label(text=f'{min_text}:.2f', x=-40, y=15, y_units='screen',
                 x_units="screen", text_align="left", text_baseline="middle", text_font_size="12px")
plot.add_layout(tick_min, "right")
# 4. max value
tick_max = Label(text=f'{max_text}:.2f', x=-50, y=510, y_units='screen',
                 x_units="screen", text_align="left", text_baseline="middle", text_font_size="12px")
plot.add_layout(tick_max, "right")

# If we plot by categorical obs, we hide all the 4 components
if Z_AXIS == "categorical facet":
    color_bar.visible = False
    colorbar_text.visible = False
    tick_min.visible = False
    tick_max.visible = False

# initial set to x axis label and y axis label
x_label = ""
y_label = ""




##spinner test##
#thanks https://discourse.bokeh.org/t/small-example-on-a-loading-div/9058/6
spinner_text = """
<!-- https://www.w3schools.com/howto/howto_css_loader.asp -->
<div class="loader">
<style scoped>
.loader {
    border: 16px solid #f3f3f3; /* Light grey */
    border-top: 16px solid #3498db; /* Blue */
    border-radius: 50%;
    width: 120px;
    height: 120px;
    animation: spin 2s linear infinite;
    displat: flex;
    margin: 250px 0 0 300px;
}

@keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
} 
</style>
</div>
"""
div_spinner = Div(text=spinner_text,visible=True)



def remove_spinner():
    div_spinner.visible = False
    plot.visible = True

def cb_update_plot(attr, old, new, type_change):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, index_cmap, X_AXIS, Y_AXIS, widget_axes, x_label, y_label, data, color_bar, colorbar_text, tick_min, tick_max
    # change type: if the button clicked needs to have an update in plotting or not
    # if type_change not in ["XYAXIS", "regression", "categorical", "geneX", "geneY", "num_facetX", "num_facetY","num_facetZ"]:
    #     curdoc().unhold()
    #     return
    
    div_spinner.visible = True
    plot.visible = False
    if type_change == "cosmetics":
        for glyph_renderer in plot.renderers:
            glyph_renderer.glyph.size = w_size_slider.value
            glyph_renderer.glyph.fill_alpha = w_alpha_slider.value
            curdoc().add_next_tick_callback(remove_spinner)

        curdoc().unhold()
        return

    data = get_data()
    #no numerical facets. => need to add which views have numerical and which dont : updates widgets

    dataset_id, dataset = get_dataset()

    facet = w_facet.value

    unique_obs = get_unique_obs(data)

    if len(unique_obs) < 3:
        palette = ["blue","red"]
    else:
        palette = Category20[len(unique_obs)]
    index_cmap = factor_cmap('obs', palette, unique_obs)

    # empty the plot glyphs and legends:
    plot.renderers = []
    if plot.legend:
        plot.legend.items = []

    # set up new x and y axis to be plotted
    X_AXIS = "geneX" if w_x_axis_radio.active == GENE_OPTION else "num_facetX"
    Y_AXIS = "geneY" if w_y_axis_radio.active == GENE_OPTION else "num_facetY"
    Z_AXIS = LABELS_GROUPING[w_category_radio.active]

    # set colorbar (the 4 components) to be invisible
    color_bar.visible = False
    colorbar_text.visible = False
    tick_min.visible = False
    tick_max.visible = False
    ALPHA =  w_alpha_slider.value
    SIZE = w_size_slider.value
    # plot again:
    if Z_AXIS == "categorical facet":
        for index, obs in enumerate(unique_obs):

            new_source = ColumnDataSource(data.loc[(data.obs == obs)])

            if len(w_regression.active) == 1:
                # plot.update(output_backend = "canvas")
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    new_source.data[X_AXIS], y=new_source.data[Y_AXIS])
                y_predicted = [round(slope, 3)*i + round(intercept, 3)
                               for i in new_source.data[X_AXIS]]

                plot.line(new_source.data[X_AXIS], np.array(y_predicted), color=index_cmap["transform"].palette[index],
                          legend_label=f"{obs}: y={str(round(slope,2))}x+{str(round(intercept,2))}, p-value={'{:e}'.format(p_value)}, r-value={str(round(r_value,2))}")
            else:
                plot.scatter(x=X_AXIS, y=Y_AXIS, source=new_source,  legend_label=obs,
                             fill_alpha=ALPHA, size=SIZE, width=0, fill_color=index_cmap["transform"].palette[index])
        plot.legend.location = "top_right"
        plot.legend.click_policy = "hide"
    elif Z_AXIS == "gene":
        mapper = get_mapper()
        plot.scatter(x=X_AXIS, y=Y_AXIS, source=ColumnDataSource(data),
                     fill_alpha=ALPHA, size=SIZE, width=0, fill_color={'field': 'geneZ', 'transform': mapper})

        # change component values texts
        colorbar_text.text = f'{w_gene3.value}'
        tick_min.text = f'{data["geneZ"].min():.2f}'
        tick_max.text = f'{data["geneZ"].max():.2f}'

        # make them visible again
        color_bar.visible = True
        colorbar_text.visible = True
        tick_min.visible = True
        tick_max.visible = True
    else:
        mapper = get_mapper()
        plot.scatter(x=X_AXIS, y=Y_AXIS, source=ColumnDataSource(data),
                     fill_alpha=ALPHA, size=SIZE, width=0, fill_color={'field': 'num_facetZ', 'transform': mapper})

        # change component values texts
        colorbar_text.text = f'{w_facet_numerical_3.value}'
        tick_min.text = f'{data["num_facetZ"].min():.2f}'
        tick_max.text = f'{data["num_facetZ"].max():.2f}'

        # make them visible again
        color_bar.visible = True
        colorbar_text.visible = True
        tick_min.visible = True
        tick_max.visible = True

    # title text change
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

    plot.title.text = (f"{x_label} vs {y_label} - "
                       f"({dataset_id}) {title}...")

    # fetch the name of the x axis label and y axis label from the widget value.
    for ind, widget in enumerate(widget_axes):
        if X_AXIS in widget.name:
            x_label = widget_axes[ind].value
        if Y_AXIS in widget.name:
            y_label = widget_axes[ind].value
    plot.xaxis.axis_label = f"{x_label}"
    plot.yaxis.axis_label = f"{y_label}"

    w_download_filename.text = f"exp_{dataset_id}_{facet}_{x_label}_{y_label}.tsv"
    curdoc().add_next_tick_callback(remove_spinner)
    # div_spinner.visible = False

    curdoc().unhold()


# convenience shortcut
update_plot = partial(cb_update_plot, attr=None, old=None,
                      new=None, type_change="XYAXIS")

# run it directly to ensure there are initial values
update_plot()


def cb_dataset_change(attr, old, new):
    """Dataset change."""
    lg.info("Dataset Change to:" + new)
    curdoc().hold()
    update_facets()
    update_genes()
    update_numerical_facets()
    update_plot()


cb_download = CustomJS(
    args=dict(data=ColumnDataSource(data).data,
              columns=[x for x in ColumnDataSource(
                  data).data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")

# not active for now. removed from root()
w_regression.on_change("active", partial(
    cb_update_plot, type_change="regression"))

# two options for coloring
w_facet.on_change("value", partial(cb_update_plot, type_change="categorical"))
w_gene3.on_change("value", partial(cb_update_plot, type_change="categorical"))
w_facet_numerical_3.on_change("value",partial(cb_update_plot, type_change="categorical"))
# widgets responsible for changing the x and y axis options
w_gene1.on_change("value", partial(cb_update_plot, type_change="geneX"))
w_gene2.on_change("value", partial(cb_update_plot, type_change="geneY"))
w_facet_numerical_1.on_change("value", partial(cb_update_plot, type_change="num_facetX"))
w_facet_numerical_2.on_change("value", partial(cb_update_plot, type_change="num_facetY"))

# widgets for changing which one of the (gene,gene,facet,facet) are we putting on x and y axis
w_x_axis_radio.on_change("active", partial(
    cb_update_plot, type_change="XYAXIS"))
w_y_axis_radio.on_change("active", partial(
    cb_update_plot, type_change="XYAXIS"))
w_category_radio.on_change("active", partial(
    cb_update_plot, type_change="categorical"))

w_alpha_slider.on_change("value", partial(cb_update_plot, type_change="cosmetics"))
w_size_slider.on_change("value", partial(cb_update_plot, type_change="cosmetics"))

w_dataset_id.on_change("value", cb_dataset_change)
w_download.js_on_click(cb_download)

cb = CustomJS(args=dict(div_spinner=div_spinner)
              ,code='''
              console.log('cb triggered!')
              div_spinner.change.emit()''')
              
div_spinner.js_on_change('visible',cb)
# curdoc().add_root(
#     column([
#         row([
#             column([
#                 row([w_gene1, w_gene2],
#                     sizing_mode='stretch_width'),
#                 row([w_facet_numerical_1, w_facet_numerical_2],
#                     sizing_mode='stretch_width'),
#                 row([w_x_axis_radio, w_y_axis_radio],
#                     sizing_mode='stretch_width'),
#                 row([w_gene3, w_facet, w_facet_numerical_3, w_category_radio,w_download],
#                     sizing_mode='stretch_width'),
#                 row([w_alpha_slider, w_size_slider],
#                     sizing_mode='stretch_width')],
#                    sizing_mode='stretch_width')],
#             sizing_mode='stretch_width'),
#         row([w_div_title_author],
#             sizing_mode='stretch_width'),
#         # TODO regression doesnt work with webgl
#         row([w_gene_not_found],
#             sizing_mode='stretch_width'),
#         row([plot],
#             sizing_mode='stretch_width', name="main_layout"),
#         row([w_dataset_id],
#             sizing_mode='stretch_width'),
#     ], sizing_mode='stretch_width')
# )



curdoc().add_root(row([
            column([
                row([w_gene1,w_gene2]),
                row([w_facet_numerical_1,w_facet_numerical_2]),
                row([w_x_axis_radio,w_y_axis_radio]),
                row([w_gene3,w_facet]),
                row([w_category_radio,w_download]),
                row([w_alpha_slider,w_size_slider]),
                row([w_div_title_author]),
                row([w_dataset_id]),
            ]),
            column([
                row([plot,div_spinner], sizing_mode="stretch_width"),
                ]),


    ])
)