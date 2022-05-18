"""Simple expression visualization."""

from functools import partial
import sqlite3

import pandas as pd

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource
from bokeh.models import DataTable, TableColumn, ScientificFormatter
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Select, TextInput, Div, Button
from bokeh.plotting import figure, curdoc

from beehive import config, util


curdoc().template_variables['config'] = config
curdoc().template_variables['view_name'] = 'Gene Expression'


create_widget = partial(util.create_widget, curdoc=curdoc())


datadir = util.get_datadir("gene_expression")
dbfile = datadir / "gene_expression.db"
assert dbfile.exists()
db = sqlite3.connect(dbfile)

args = curdoc().session_context.request.arguments

datasets = pd.read_sql('SELECT * FROM STUDY', db)
datasets['full'] = datasets['title'] + ', ' + datasets['author']


# WIDGETS
def_dataset_id = util.getarg(args, 'dataset_id',
                             datasets['dataset_id'].iloc[0])
dataset_options = [(x['dataset_id'], x['full'])
                   for _, x in datasets.iterrows()]


w_div_title_author = Div(text="")

w_dataset_id = create_widget("dataset_id", Select, title="Dataset",
                             options=dataset_options,
                             default=def_dataset_id)

w_gene = create_widget("gene", TextInput, default='APOE')
w_facet = create_widget("facet", Select, options=[], title="Group by")
w_plottype = create_widget("plottype", Select, title="Show",
                           options=["boxplot", "mean/std"])

w_download = Button(label='Download', align='end')
w_download_filename = Div(text="", visible=Fals,
                          name="download_filename")

# To display text if the gene is not found
w_gene_not_found = Div(text="")


#
# Data handling & updating interface
#


def get_facets():
    """Get available facets for a dataset."""
    dataset_id = w_dataset_id.value
    sql = f'''SELECT DISTINCT(cat_name)
                FROM data
               WHERE dataset_id == "{dataset_id}"'''
    facets = list(sorted(pd.read_sql(sql, db)['cat_name']))
    return facets


def get_data() -> pd.DataFrame:
    """Retrieve data from a dataset, gene & facet."""
    dataset_id = w_dataset_id.value
    gene = w_gene.value
    facet = w_facet.value
    sql = f'''SELECT *
                FROM data
               WHERE dataset_id == "{dataset_id}"
                 AND cat_name == "{facet}"
                 AND gene == "{gene}"
           '''
    data = pd.read_sql(sql, db)
    data.sort_index(axis=1, inplace=True)

    # default settings for a boxplot -
    # override if other further down
    data['_segment_top'] = data['q99']
    data['_bar_top'] = data['q80']
    data['_bar_median'] = data['q50']
    data['_bar_bottom'] = data['q10']
    data['_segment_bottom'] = data['q01']
    # print(data.head(2).T)
    return data


def get_categories() -> list:
    """Return possible categories for a facet."""
    dataset_id = w_dataset_id.value
    facet = w_facet.value
    sql = f'''SELECT DISTINCT("cat_value")
                FROM data
               WHERE dataset_id == "{dataset_id}"
                 AND cat_name == "{facet}"
           '''
    categs = list(sorted(pd.read_sql(sql, db)['cat_value']))
    return categs


def get_dataset() -> pd.Series:
    """Return the current dataset record."""
    dataset_id = w_dataset_id.value
    rv = datasets[datasets['dataset_id'] == dataset_id]
    assert len(rv) == 1
    return rv.iloc[0]


#
# Change & Initialize interface
#
def update_facets():
    """Update interface for a specific dataset."""
    facets = get_facets()
    w_facet.options = facets
    if w_facet.value not in facets:
        w_facet.value = facets[0]


update_facets()


#
# Create plot
#
plot = figure(background_fill_color="#efefef", x_range=[],
              plot_height=400, title="Plot",
              toolbar_location='right')

source = ColumnDataSource(get_data())

table = DataTable(source=source,
                  margin=10,
                  index_position=None,
                  columns=[
                      TableColumn(field='cat_value', title='Category'),
                      TableColumn(field='mean', title='Mean',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q50', title='Median',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q01', title='1% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q20', title='20% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q80', title='80% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                      TableColumn(field='q99', title='99% Quantile',
                                  formatter=ScientificFormatter(precision=2)),
                  ])


# create segments
elements = dict(
    vbar=plot.vbar(source=source, x='cat_value', top='_bar_top',
                   bottom='_bar_bottom', width=0.85, name="barplot",
                   fill_color=config.colors.color6,
                   line_color="black"),
    seg_v_up=plot.segment(source=source, x0='cat_value', x1='cat_value',
                          y0='_bar_top', y1='_segment_top', line_color='black'),
    seg_h_up=plot.rect(source=source, x='cat_value', height=0.001,
                       y='_segment_top', width=0.4, line_color='black'),
    seg_v_dw=plot.segment(source=source, x0='cat_value', x1='cat_value',
                          y0='_segment_bottom', y1='_bar_bottom',
                          line_color='black'),
    seg_h_dw=plot.rect(source=source, x='cat_value', height=0.001,
                       y='_segment_bottom', width=0.4, line_color='black'),
    seg_h_med=plot.rect(source=source, x='cat_value', height=0.001,
                        y='_bar_median', width=0.85, line_width=2,
                        line_color='black'),
)


def update_plot():
    """Update the plot."""
    global plot, source
    data = get_data()
    dataset = get_dataset()
    facet = w_facet.value
    gene = w_gene.value

    w_div_title_author.text = \
        f"""
        <dl>
          <dt>Title:</dt><dd>{dataset['title']}</dd>
          <dt>Author:</dt><dd>{dataset['author']}</dd>
        </dl>
        """
    gene = w_gene.value
    plot.x_range.factors = list(data['cat_value'])
    w_download_filename.text = f"exp_{dataset['dataset_id']}_{facet}_{gene}.tsv"

    if w_plottype.value == "boxplot":
        pttext = 'Boxplot'
    else:
        data['_segment_top'] = data['mean'] + data['std']
        data['_bar_top'] = data['mean']
        data['_bar_bottom'] = 0
        data['_segment_bottom'] = 0
        data['_bar_median'] = data['mean']
        pttext = 'Mean/std'

    ymax = max(1, 1.01 * data['_segment_top'].max())
    source.data = data

    title = dataset['title'][:60]
    plot.title.text = (f"{pttext} {gene} vs {facet} - "
                       f"({dataset.dataset_id}) {title}...")
    plot.y_range.end = ymax


update_plot()


#
# widget callbacks
#
def _dataset_change(attr, old, new):
    """Dataaset change."""
    update_facets()
    update_plot()


def _gene_change(attr, old, new):
    """Gene change callback."""
    # First check if the gene is in the database
    dataset_id = w_dataset_id.value
    gene = w_gene.value
    sql = f"""SELECT *
                FROM data
               WHERE dataset_id = "{dataset_id}"
                 AND gene = "{gene}"
               LIMIT 1 """

    check_gene_in_db = pd.read_sql(sql, db)
    if len(check_gene_in_db) == 0:

        # hide the plot and show the help div
        plot.visible = False
        table.visible = False
        w_gene_not_found.visible = True

        # find candidates
        gene_find = gene
        while len(gene_find) > 2:
            sql = \
                f"""
                SELECT DISTINCT gene FROM data
                 WHERE dataset_id = "{dataset_id}"
                   AND gene LIKE "%{gene_find}%"
                 LIMIT 50
                """
            genes = list(pd.read_sql(sql, db)['gene'])

            def atag(gene):
                """Inject javascript call to update gene widget."""
                return (f'''
                    <a onclick='Bokeh.documents[0].get_model_by_name("gene").value = "{gene}";'>
                       {gene}
                    </a>''')  # noqa: E501

            if len(genes) > 0:
                genestxt = "<p>Did you mean: " + \
                    ", ".join(map(atag, genes)) + "?"
                import textwrap
                genestxt = "\n".join(textwrap.wrap(genestxt))
                break
            else:
                gene_find = gene_find[:-1]
        else:
            # Have not found candidate?
            genestxt = "\nCannot find a candidate."

        w_gene_not_found.text = (
            f"<b>Gene not found: {gene}!</b>{genestxt}")
    else:
        plot.visible = True
        table.visible = True
        w_gene_not_found.visible = False
        update_plot()


def _update_plot(attr, old, new):
    update_plot()


download_callback = CustomJS(
    args=dict(data=source.data,
              columns=[x for x in source.data.keys() if not x.startswith('_')],
              filename_div=w_download_filename),
    code="exportToTsv(data, columns, filename_div.text);")

w_gene.on_change("value", _gene_change)
w_dataset_id.on_change("value", _dataset_change)
w_facet.on_change("value", _update_plot)
w_plottype.on_change("value", _update_plot)
w_download.js_on_click(download_callback)


#
# Build the document
#
curdoc().add_root(
    column([
        row([w_dataset_id, w_download_filename],
            sizing_mode='stretch_width'),
        row([w_gene, w_facet, w_plottype, w_download],
            sizing_mode='stretch_width'),
        row([w_div_title_author],
            sizing_mode='stretch_width'),
        row([w_gene_not_found],
            sizing_mode='stretch_width'),
        row([plot],
            sizing_mode='stretch_width'),
        row([table],
            sizing_mode='stretch_width'),
    ], sizing_mode='stretch_width')
)
