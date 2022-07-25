import logging
from functools import partial
import logging
from pprint import pprint
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource,CheckboxGroup,Band
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, TextInput, Div,
                                  Button, AutocompleteInput)
from bokeh.plotting import figure, curdoc
from bokeh.transform import factor_cmap
from bokeh.palettes import Category20

from beehive import config, util, expset

lg = logging.getLogger('ScatterExpression')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Scatter Expression'


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


## note: human/plant has meta as 'models' , mouse has meta as 'model'
siblings = expset.get_dataset_siblings(w_dataset_id.value)

#making widgets for the following:

#gene, has an autocomplete feature.. write in text
#4 widgets, determine what to show?
#gene vs gene
#numerical facet vs numerical facet
#numerical facet vs gene etc...
w_gene1 = create_widget("geneX", AutocompleteInput,
                       completions=[], default='APOE', title="Gene X")
w_gene2 = create_widget("geneY", AutocompleteInput,
                       completions=[], default='TREM2', title = "Gene Y")
w_facet_numerical_1 = create_widget("num_facetX",Select, 
                    options=[], title="Select Numerical Facet on X")
w_facet_numerical_2 = create_widget("num_facetY",Select, 
                        options=[], title="Select Numerical Facet on Y")

widget_axes = [w_gene1,w_gene2,w_facet_numerical_1,w_facet_numerical_2]

w_regression  = CheckboxGroup(labels=["Regression Lines"], active=[])

FIXED_OPTIONS = [("geneX","Gene X"),("geneY","Gene Y"),("num_facetX","Numerical Facet on X"),("num_facetY","Numerical Facet on Y")]

w_x_axis = create_widget("x_axis",Select, 
                        options=FIXED_OPTIONS, title="Select X-Axis",default = 'geneX')
w_y_axis =  create_widget("y_axis",Select, 
                        options=FIXED_OPTIONS, title="Select Y-Axis",default = 'geneY')

#categorical facet grouping of data...
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
    ##TODO maybe check for numerical/categorical here?
    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene1}")

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene2}")

    data = pd.DataFrame(dict(
        geneX = expset.get_gene(dataset_id, gene1)[:,0],
        geneY = expset.get_gene(dataset_id, gene2)[:,0],
        obs = expset.get_meta(dataset_id, facet)[:,0],
        num_facetX = expset.get_meta(dataset_id,num_facet1,raw=True)[:,0],
        num_facetY = expset.get_meta(dataset_id,num_facet2,raw=True)[:,0]
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

plot = figure()
data = get_data()

source = ColumnDataSource(data)
unique_obs = get_unique_obs(data)
dataset_id, dataset = get_dataset()


index_cmap = factor_cmap('obs', Category20[len(unique_obs)], unique_obs)


glyphs = []
sources = []

X_AXIS = "geneX"
Y_AXIS = "geneY"

##need to have multiple glyphs not a single multi glyph:
for index,obs in enumerate(unique_obs):

    sourcedf = pd.DataFrame(source.data)
    new_source = ColumnDataSource(sourcedf.loc[(sourcedf.obs == obs)])
    sources = sources + [new_source]
    glyph = plot.scatter(x=X_AXIS, y=Y_AXIS, source=new_source,  legend_label=obs,
    fill_alpha=0.7, size=5,width=0, fill_color = index_cmap["transform"].palette[index])

    glyphs = glyphs + [glyph]

plot.legend.location = "top_right"
plot.legend.click_policy = "hide"
x_label = ""
y_label = ""


# def cb_update_plot(attr, old, new,type_change,axis):
def cb_update_plot(attr, old, new, type_change):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, sources, index_cmap,glyphs,source, X_AXIS, Y_AXIS, widget_axes,x_label,y_label
    if type_change != w_x_axis.value and type_change != w_y_axis.value and type_change not in ["XYAXIS","regression","categorical"]:
        curdoc().unhold()
        return

    data = get_data()
    dataset_id, dataset = get_dataset()

    facet = w_facet.value

    source = ColumnDataSource(data)


    unique_obs = get_unique_obs(data)

    #TODO fix None => strings
    index_cmap = factor_cmap('obs', Category20[len(unique_obs)], unique_obs)
    
    plot.renderers = []
    plot.legend.items = []

    glyphs = []
    sources = []

    X_AXIS = w_x_axis.value
    Y_AXIS = w_y_axis.value
    for index,obs in enumerate(unique_obs):

        sourcedf = pd.DataFrame(source.data)
        new_source = ColumnDataSource(sourcedf.loc[(sourcedf.obs == obs)])
        sources = sources + [new_source]

        if len(w_regression.active) == 1:

            slope, intercept, r_value, p_value, std_err = stats.linregress(new_source.data[X_AXIS],y = new_source.data[Y_AXIS])
            y_predicted = [slope*i + intercept  for i in new_source.data[X_AXIS]]
            glyph = plot.line(new_source.data[X_AXIS],y_predicted,color=index_cmap["transform"].palette[index],
            legend_label=f"{obs}: y={str(round(slope,2))}x+{str(round(intercept,2))}, p-value={'{:e}'.format(p_value)}, r-value={str(round(r_value,2))}")

            ##TRYING CONF Intervals
            # sem = lambda x: x.std() / np.sqrt(x.size)
            # df2 = pd.DataFrame(new_source.data)[Y_AXIS].rolling(window=100).agg({"y_mean": np.mean, "y_std": np.std, "y_sem": sem})

            # df2 = pd.concat([pd.DataFrame(new_source.data)[X_AXIS], df2], axis=1)
            # df2 = df2.fillna(method='bfill')
            # df2['lower'] = df2['y_mean'] - df2['y_std']
            # df2['upper'] = df2['y_mean'] + df2['y_std']

            # print(df2)
            # band = Band(base = X_AXIS,lower="lower", 
            # upper="upper", source= ColumnDataSource(df2), level='underlay',fill_alpha=1.0, line_width=1, 
            # line_color=index_cmap["transform"].palette[index])
            # if index == 1:
            #     plot.add_layout(band)

            # band = plot.line(x = X_AXIS, y = "lower" ,source= ColumnDataSource(df2), line_width=1, 
            # line_color=index_cmap["transform"].palette[index])


        else:
            glyph = plot.scatter(x=X_AXIS, y=Y_AXIS, source=new_source,  legend_label=obs,
            fill_alpha=0.7, size=5,width=0, fill_color = index_cmap["transform"].palette[index])

        glyphs = glyphs + [glyph]

    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"

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
    update_plot()


def cb_sibling_change(attr, old, new):
    lg.debug("Sibling change: " + new)
    w_dataset_id.value = new

##TODO: fix download button. only downloads old data
def cb_download():
    data = get_data()
    source.data = data
    return CustomJS(
    args=dict(data=source.data,
              columns=[x for x in [X_AXIS,Y_AXIS,"obs"] if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")



w_regression.on_change("active",partial(cb_update_plot,type_change = "regression"))
w_x_axis.on_change("value",partial(cb_update_plot,type_change = "XYAXIS"))
w_y_axis.on_change("value",partial(cb_update_plot,type_change = "XYAXIS"))
w_facet.on_change("value", partial(cb_update_plot,type_change="categorical"))

w_gene1.on_change("value",partial(cb_update_plot,type_change = "geneX"))
w_gene2.on_change("value", partial(cb_update_plot,type_change="geneY"))
w_facet_numerical_1.on_change("value", partial(cb_update_plot,type_change="num_facetX"))
w_facet_numerical_2.on_change("value", partial(cb_update_plot,type_change="num_facetY"))


w_sibling.on_change("value", cb_sibling_change)
w_dataset_id.on_change("value", cb_dataset_change)
w_download.js_on_click(cb_download())


curdoc().add_root(
    column([
        row([
            column([
                row([w_gene1,w_gene2,w_facet], 
                    sizing_mode='stretch_width'),
                row([w_facet_numerical_1,w_facet_numerical_2,w_sibling],
                    sizing_mode='stretch_width'),
                row([w_x_axis,w_y_axis,w_download],
                    sizing_mode='stretch_width')],   
                sizing_mode='stretch_width')],
            sizing_mode='stretch_width'),
        row([w_div_title_author],
            sizing_mode='stretch_width'),
        row([w_regression]),
        row([w_gene_not_found],
            sizing_mode='stretch_width'),
        row([plot],
            sizing_mode='stretch_width'),
        row([w_dataset_id, w_download_filename],),
    ], sizing_mode='stretch_width')
)