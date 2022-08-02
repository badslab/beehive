import logging
from functools import partial
import logging
from venv import create
import pandas as pd
import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, RadioGroup, CheckboxGroup, Slider,MultiSelect, LinearColorMapper
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Select, Div,
                                  Button, AutocompleteInput)
from bokeh.plotting import figure, curdoc
from bokeh.models.transforms import CustomJSTransform


from beehive import config, util, expset
from bokeh.util.hex import cartesian_to_axial

lg = logging.getLogger('ScatterExpression')
lg.setLevel(logging.DEBUG)
lg.info("startup")

curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Hexbin Expression'


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
LABELS_AXIS = ["Gene", "Numerical Facet"]
w_gene1 = create_widget("gene1", AutocompleteInput,
                       completions=[], default='APOE', title = "Gene")
w_gene2 = create_widget("gene2", AutocompleteInput,
                       completions=[], default='TREM2', title = "Gene")
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
w_gene3 = create_widget("gene3", AutocompleteInput,
                       completions=[], default='TREM2', title = "Gene")
#color by continuous facet
w_facet_numerical_3 = create_widget("num_facetZ", Select,
                                    options=[], title="Numerical Facet")
#color/ yes or no?
w_z_axis = CheckboxGroup(labels=["Color by Z-Axis"], active=[0])
#color by gene or numerical facet?
w_category_radio = create_widget(
    "category_radio", RadioGroup, labels=LABELS_GROUPING, default=0, title="Grouped:", value_type=int)

#user selecting alpha
w_alpha_slider = create_widget("alpha_picker",Slider,start=0.2, end=1, default=0.5,step=0.05,title = "Opacity",value_type = float)

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
    facet = w_facet.value
    ##TODO maybe check for numerical/categorical here?
    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene1}")

    lg.warning(f"!! Getting data for {dataset_id} {facet} {gene2}")

    ##TODO scheme to get data => based on what's selected only
    data = pd.DataFrame(dict(
        gene1 = expset.get_gene(dataset_id, gene1)[:,0],
        gene2 = expset.get_gene(dataset_id, gene2)[:,0],
        gene3 = expset.get_gene(dataset_id, gene3)[:,0],
        obs=expset.get_meta(dataset_id, facet)[:, 0],
        ))

    return data

def get_dataset():
    """Return the current dataset id and record."""
    dataset_id = w_dataset_id.value
    return dataset_id, datasets[dataset_id]

def get_unique_obs(data):
    """Fetch the unique facets from the data"""
    unique_obs = pd.DataFrame(data)['obs'].unique()
    w_subset_select.options = unique_obs.tolist()

#
# Create plot
#
def modify_data():
    data = get_data()
    global unique_obs
    unique_obs = get_unique_obs(data)
    q, r = cartesian_to_axial(data["gene1"], data["gene2"], 0.1, "pointytop")
    df = pd.DataFrame(dict(r=r,q=q,avg_exp = None))
    groups = df.groupby(["q","r"])
    dicts = []
    ###assign axial coordinates to get average gene expression###
    for key,val in groups.groups.items():
        q,r = key
        counts = len(val)
        avg_exp = data['gene3'].loc[groups.groups.get((q,r))].mean()
        dicts = dicts + [{"q" : q,"r" : r,"counts" : counts, "avg_exp" :avg_exp,"index_list":val}]

    final_result = pd.DataFrame(dicts)

    return final_result


plot = figure()
dataset_id, dataset = get_dataset()
data = modify_data()


#TODO make 1 and 99 percentiles
mapper = LinearColorMapper(
    palette='Magma256',
    low=data["avg_exp"].max(),
    high=data["avg_exp"].min())

v_func_alpha  = """
var new_xs = new Array(xs.length)
for(var i = 0; i < xs.length; i++) {
    new_xs[i] = alpha_map[xs[i]]
}
return new_xs
"""

source = ColumnDataSource(data)

length = len(sorted(np.unique(source.data["counts"]).tolist()))
alpha_map = dict(zip(sorted(np.unique(source.data["counts"]).tolist()),
[w_alpha_slider.value if x/length < w_alpha_slider.value else x/length for x in range(0,length)]))

numerical_alpha_transformer = CustomJSTransform(args={"alpha_map": alpha_map}, v_func=v_func_alpha)

#TODO: what x and y are we plotting

#TODO: do we need to color? (add another hextile or not.)

#TODO: if yes, what color sceheme are we using? facet or gene expression

#TODO: show subsets only => check if subsets are selected, if not, then show all.

#TODO: extras: add colorbar, add xy title labels, legend

X_AXIS = "gene1"
Y_AXIS = "gene2"


x_label = ""
y_label = ""

def cb_update_plot(attr, old, new,type_change,axis):
    """Populate and update the plot."""
    curdoc().hold()
    global plot, index_cmap,source, X_AXIS, Y_AXIS, widget_axes,x_label,y_label,merged_image
    data = modify_data()
    dataset_id, dataset = get_dataset()
    facet = w_facet.value
    source = ColumnDataSource(data)

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
update_plot = partial(cb_update_plot, attr=None, old=None, new=None,type_change=None,axis=None)

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

w_gene1.on_change("value",partial(cb_update_plot,type_change = "gene1",axis="x"))
w_gene2.on_change("value", partial(cb_update_plot,type_change="gene2",axis="y"))
w_facet.on_change("value", partial(cb_update_plot,type_change=None,axis=None))
w_facet_numerical_1.on_change("value", partial(cb_update_plot,type_change="num_facet1",axis="x"))
w_facet_numerical_2.on_change("value", partial(cb_update_plot,type_change="num_facet2",axis="y"))

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
                row([w_x_axis_radio,w_y_axis_radio,w_category_radio],
                    sizing_mode='stretch_width'),
                row([w_facet,w_subset_select, w_alpha_slider,w_download],
                    sizing_mode='stretch_both'),
                row([w_z_axis],
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