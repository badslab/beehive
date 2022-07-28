import logging
from functools import partial
import logging
import pandas as pd
from scipy import stats
import numpy as np
from bokeh.layouts import column, row,layout
from bokeh.models import ColumnDataSource,CheckboxGroup, RadioGroup, Slope, LinearColorMapper, ColorBar
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


#TODO
#clear cookies for user/password? ==> in auth not here..
#print(curdoc().session_context.request.cookies)

w_div_title_author = Div(text="")

#datasets with titles

# Dataset
dataset_options = [(k, "{short_title}, {short_author}".format(**v))
                   for k, v in datasets.items()]

##TODO setting manually
DATASET_NUMBER = 3

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=dataset_options[DATASET_NUMBER][0],
                             visible=True,)

# Possible siblings of this dataset
# siblings = expset.get_dataset_siblings(w_dataset_id.value)
# sibling_options = []
# for k, v in siblings.items():
#     sname = f"{v['organism']} / {v['datatype']}"
#     sibling_options.append((k, sname))

# w_sibling = create_widget("view", Select,
#                           options=sibling_options,
#                           default=w_dataset_id.value,
#                           update_url=False)


## note: human/plant has meta as 'models' , mouse has meta as 'model'
# siblings = expset.get_dataset_siblings(w_dataset_id.value)

#making widgets for the following:

#gene, has an autocomplete feature.. write in text
#4 widgets, determine what to show?
#gene vs gene
#numerical facet vs numerical facet
#numerical facet vs gene etc...
w_gene1 = create_widget("geneX", AutocompleteInput,
                       completions=[], default='APOE', title="Gene", case_sensitive = False)
w_gene2 = create_widget("geneY", AutocompleteInput,
                       completions=[], default='TREM2', title = "Gene", case_sensitive = False)
w_facet_numerical_1 = create_widget("num_facetX",Select, 
                    options=[], title="Numerical Facet")
w_facet_numerical_2 = create_widget("num_facetY",Select, 
                        options=[], title="Numerical Facet")

widget_axes = [w_gene1,w_gene2,w_facet_numerical_1,w_facet_numerical_2]

w_regression  = CheckboxGroup(labels=["Regression Lines"], active=[])

FIXED_OPTIONS = [("geneX","Gene X"),("geneY","Gene Y"),("num_facetX","Numerical Facet on X"),("num_facetY","Numerical Facet on Y")]

# w_x_axis = create_widget("x_axis",Select, 
#                         options=FIXED_OPTIONS, title="Select X-Axis",default = 'geneX')
# w_y_axis =  create_widget("y_axis",Select, 
#                         options=FIXED_OPTIONS, title="Select Y-Axis",default = 'geneY')

LABELS_AXIS = ["Gene","Numerical Facet"]
LABELS_GROUPING = ["facet","gene"]

w_x_axis_radio = create_widget("x_axis_radio",RadioGroup,labels = LABELS_AXIS, default = 0, title = "X-Axis",value_type = int)
w_y_axis_radio = create_widget("y_axis_radio",RadioGroup,labels = LABELS_AXIS, default = 0, title = "Y-Axis", value_type = int)

#categorical facet grouping of data... or gene3 grouping
w_facet = create_widget("facet", Select, options=[], title="Group by Categories")

w_gene3 = create_widget("gene_categorical", AutocompleteInput,
                       completions=[], default="APOE", title="Group by Gene Expression", case_sensitive = False)

w_category_radio = create_widget("category_radio",RadioGroup,labels = LABELS_GROUPING, default = 0, title = "Grouped:",value_type = int)

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
    options = expset.get_facet_options_numerical(w_dataset_id.value)
    w_facet_numerical_1.options = options
    w_facet_numerical_2.options = options
    if w_facet_numerical_1.value not in [x[0] for x in options]:
        # set a default
        w_facet_numerical_1.value = options[0][0]
    if w_facet_numerical_2.value not in [x[0] for x in options]:
        # set a default
        w_facet_numerical_2.value = options[0][0]

update_numerical_facets()


def get_data() -> pd.DataFrame:
    """Retrieve data from a dataset, gene & facet."""
    dataset_id = w_dataset_id.value
    gene1 = w_gene1.value
    gene2 = w_gene2.value
    facet = w_facet.value
    num_facet1 = w_facet_numerical_1.value
    num_facet2 = w_facet_numerical_2.value
    # x_axis = w_x_axis.value
    # y_axis = w_y_axis.value

    geneX = None
    geneY = None
    num_facetX = None
    num_facetY = None
    geneZ = None
    x_axis = w_x_axis_radio.active
    y_axis = w_y_axis_radio.active
    gene3 = w_gene3.value

    if x_axis == GENE_OPTION:
        geneX = expset.get_gene(dataset_id, gene1)[:,0]
    else:
        num_facetX = expset.get_meta(dataset_id,num_facet1,raw=True)[:,0]

    if y_axis == GENE_OPTION:
        geneY = expset.get_gene(dataset_id, gene2)[:,0]
    else:
        num_facetY = expset.get_meta(dataset_id,num_facet2,raw=True)[:,0]

    if gene3:
        geneZ = expset.get_gene(dataset_id, gene3)[:,0]

    ##TODO maybe check for numerical/categorical here?
    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene1}")

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene2}")

    data = pd.DataFrame(dict(
        geneX = geneX,
        geneY = geneY,
        obs = expset.get_meta(dataset_id, facet)[:,0],
        num_facetX = num_facetX,
        num_facetY = num_facetY,
        geneZ = geneZ
        ))
    return data




def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]

def get_unique_obs(data):
    unique_obs = pd.DataFrame(data)['obs'].unique()
    return unique_obs


#
# Create plot
#

def get_mapper():
    mapper = LinearColorMapper(
    palette='Magma256',
    low=data["geneZ"].min(),
    high=data["geneZ"].max())
    return mapper

plot = figure(output_backend = "webgl")
data = get_data()


unique_obs = get_unique_obs(data)
dataset_id, dataset = get_dataset()


index_cmap = factor_cmap('obs', Category20[len(unique_obs)], unique_obs)



X_AXIS = "geneX"
Y_AXIS = "geneY"
Z_AXIS = "facet" if w_category_radio.active == 0 else "gene"

##need to have multiple glyphs not a single multi glyph:
if Z_AXIS == "facet":
    for index,obs in enumerate(unique_obs):

        new_source = ColumnDataSource(data.loc[(data.obs == obs)])

        plot.scatter(x=X_AXIS, y=Y_AXIS, source=new_source,  legend_label=obs,
        fill_alpha=0.7, size=5,width=0, fill_color = index_cmap["transform"].palette[index])
    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"
else:
    mapper = get_mapper()
    plot.scatter(x=X_AXIS, y=Y_AXIS, source=ColumnDataSource(data),
            fill_alpha=0.7, size=5,width=0, fill_color = {'field': 'geneZ', 'transform': mapper})
    # color_bar = ColorBar(color_mapper=mapper, location=(0,0),label_standoff=12,title=w_gene3.value)
    # plot.add_layout(color_bar,"right")





x_label = ""
y_label = ""


def cb_update_plot(attr, old, new, type_change):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, index_cmap, X_AXIS, Y_AXIS, widget_axes,x_label,y_label, data
    if type_change not in ["XYAXIS","regression","categorical","geneX","geneY","num_facetX","num_facetY"]:
        curdoc().unhold()
        return

    data = get_data()
    dataset_id, dataset = get_dataset()

    facet = w_facet.value
    
    unique_obs = get_unique_obs(data)

    #TODO fix None => strings
    index_cmap = factor_cmap('obs', Category20[len(unique_obs)], unique_obs)

    if plot.right:
        plot.right = []

    plot.renderers = []
    if plot.legend:
        plot.legend.items = []

    X_AXIS = "geneX" if w_x_axis_radio.active == GENE_OPTION else "num_facetX"
    Y_AXIS = "geneY" if w_y_axis_radio.active == GENE_OPTION else "num_facetY"
    Z_AXIS = "facet" if w_category_radio.active == 0 else "gene"

    if Z_AXIS == "facet":
        for index,obs in enumerate(unique_obs):

            new_source = ColumnDataSource(data.loc[(data.obs == obs)])

            if len(w_regression.active) == 1:
                # plot.update(output_backend = "canvas")
                slope, intercept, r_value, p_value, std_err = stats.linregress(new_source.data[X_AXIS],y = new_source.data[Y_AXIS])
                y_predicted = [round(slope,3)*i + round(intercept,3)  for i in new_source.data[X_AXIS]]
            
                # regression_line = Slope(gradient=slope, y_intercept=intercept, line_color=index_cmap["transform"].palette[index])
                # plot.add_layout(regression_line)
                plot.line(new_source.data[X_AXIS],np.array(y_predicted),color=index_cmap["transform"].palette[index],
                legend_label=f"{obs}: y={str(round(slope,2))}x+{str(round(intercept,2))}, p-value={'{:e}'.format(p_value)}, r-value={str(round(r_value,2))}")
            else:
                # plot.update(output_backend = "webgl")
                plot.scatter(x=X_AXIS, y=Y_AXIS, source=new_source,  legend_label=obs,
                fill_alpha=0.7, size=5,width=0, fill_color = index_cmap["transform"].palette[index])
   
        plot.legend.location = "top_right"
        plot.legend.click_policy = "hide"
    else:
        mapper = get_mapper()
        plot.scatter(x=X_AXIS, y=Y_AXIS, source=ColumnDataSource(data),
                fill_alpha=0.7, size=5,width=0, fill_color = {'field': 'geneZ', 'transform': mapper})
       # plot.rect(x=0.5, y='geneZ', fill_color='color', width=1, height=1, source=source)

        # color_bar = ColorBar(color_mapper=mapper, location=(0,0),label_standoff=12,title=w_gene3.value)
        # plot.add_layout(color_bar,"right")


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

    plot.xaxis.axis_label = f"{x_label}"
    plot.yaxis.axis_label = f"{y_label}"
    w_download_filename.text = f"exp_{dataset_id}_{facet}_{x_label}_{y_label}.tsv"
    


    curdoc().unhold()


# convenience shortcut
update_plot = partial(cb_update_plot, attr=None, old=None, new=None,type_change = "XYAXIS")

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
              columns=[x for x in ColumnDataSource(data).data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")


w_regression.on_change("active",partial(cb_update_plot,type_change = "regression"))
w_facet.on_change("value", partial(cb_update_plot,type_change="categorical"))
w_gene3.on_change("value",partial(cb_update_plot,type_change="categorical"))


w_gene1.on_change("value",partial(cb_update_plot,type_change = "geneX"))
w_gene2.on_change("value", partial(cb_update_plot,type_change="geneY"))
w_facet_numerical_1.on_change("value", partial(cb_update_plot,type_change="num_facetX"))
w_facet_numerical_2.on_change("value", partial(cb_update_plot,type_change="num_facetY"))

w_x_axis_radio.on_change("active",partial(cb_update_plot,type_change="XYAXIS"))
w_y_axis_radio.on_change("active",partial(cb_update_plot,type_change="XYAXIS"))
w_category_radio.on_change("active",partial(cb_update_plot,type_change="categorical"))

# w_sibling.on_change("value", cb_sibling_change)
w_dataset_id.on_change("value", cb_dataset_change)

w_download.js_on_click(cb_download)


curdoc().add_root(
    column([
        row([
            column([
                row([w_gene1,w_gene2], 
                    sizing_mode='stretch_width'),
                row([w_facet_numerical_1,w_facet_numerical_2],
                    sizing_mode='stretch_width'),
                row([w_x_axis_radio,w_y_axis_radio],
                    sizing_mode='stretch_width'),
                row([w_gene3,w_facet,w_download,w_category_radio],
                    sizing_mode='stretch_width')],   
                sizing_mode='stretch_width')],
            sizing_mode='stretch_width'),
        row([w_div_title_author],
            sizing_mode='stretch_width'),
        #TODO regression doesnt work with webgl
        row([w_gene_not_found],
            sizing_mode='stretch_width'),
        row([plot],
            sizing_mode='stretch_width'),
        row([w_dataset_id],
            sizing_mode='stretch_width'),          
    ], sizing_mode='stretch_width')
)