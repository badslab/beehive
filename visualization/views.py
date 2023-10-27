"""Visualization & plot tools."""

from copy import copy
from dataclasses import dataclass, field
import re
from typing import Optional, List, Any


import altair as alt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from termite.vis import util

# eventually speed up altair transformations?
# alt.data_transformers.enable("vegafusion")


def try_num_sort(lst: np.ndarray[Any] | List[Any]) -> List[Any]:
    """
    Attempt a numerical sort.

    Remove constant string prefixes if found.
    """
    lst = list(lst)
    found_not_str = list(filter(lambda x: not isinstance(x, str), lst))
    prefix = None

    if len(lst) > 1 and len(found_not_str) == 0:
        # only strings, more than 1
        l2 = list(sorted(lst))
        prefix = ''
        s1, s2 = l2[0], l2[-1]
        for i in range(min(len(s1), len(s2))):
            if s1[:i] == s2[:i]:
                prefix = s1[:i]
            else:
                break

    # TODO: investigate - sometimes this function might yield both numbers & strings
    #       which crash the sort?
    def _tn(x: Any) -> Any:
        """Attempt to convert to something sortable."""
        if prefix is not None:
            x = x[len(prefix):]
        if x == 'nan':
            return 'nan'
        try:
            return float(x)
        except ValueError:
            return x
    return list(sorted(lst, key=_tn))


@dataclass
class TermiteView:
    """Global dataclass for all termite views."""

    data: pd.DataFrame
    name: str = "Termite plot"
    key: Optional[str] = None
    title: Optional[str] = None
    subject: Optional[str] = None

    ctxPlot: DeltaGenerator = field(default_factory=st.empty)
    ctxData: DeltaGenerator = field(default_factory=st.container)
    ctxParam: DeltaGenerator = field(default_factory=st.container)

    def get_key(self) -> str:
        """Return of create a nice url name for this view."""
        if self.key is not None:
            return self.key
        return re.sub(r'[\(\) ]', '', self.name)

    def __post_init__(self) -> None:
        """Prepare for streamlit & check some data."""
        if self.ctxPlot is None:
            self.ctxPlot = st.empty()
        if self.ctxData is None:
            self.ctxData = st.container()

        self.ctxParam1, self.ctxParam2 = \
            self.ctxParam.columns(2)

        if self.title is None:

            title = self.name
            if self.subject is not None:
                title += f" of {self.subject}"
            self.title = title
        if isinstance(self.data, pd.Series):
            self.data = pd.DataFrame(self.data)

        self.numcol = [x for x in self.data.columns
                       if 'O' not in self.data[x].dtype.kind]

        self.catcol = [x for x in self.data.columns
                       if 'O' in self.data[x].dtype.kind]

    def get_catcol(self,
                   colname: str = '1') -> str:
        """Return one categorical Column."""
        if len(self.catcol) == 1:
            return self.catcol[0]
        else:
            return util.selectbox_mem(
                label=f"Categorical Column {colname}",
                key=f'view_catcol_{colname}',
                options=self.catcol,
                context=self.ctxParam1)

    def get_numcol(self,
                   colname: str = '1',
                   exclude: List[str] = []) -> str:
        """Return name of a numeric column (and gene!)."""
        options = copy(self.numcol)
        for ex in exclude:
            while ex in options:
                options.remove(ex)

        if len(self.numcol) == 1:
            return self.numcol[0]
        else:
            return util.selectbox_mem(
                label=f"Numerical Column {colname}",
                key=f'view_numcol_{colname}',
                options=options,
                context=self.ctxParam2,)


class Multi():
    """Collect mulitple view, show selectbox, show one."""

    def __init__(self,
                 *plots: TermiteView,
                 context:  DeltaGenerator = st.sidebar):
        """Create multiview.

        * Store context to render streamlit object.
        * prepare plots & store names for the select objects
          later.
        """
        self.context = context
        self.plots = {a.get_key(): a for a in plots}
        self.plotnames = {
            a.get_key(): a.name
            for a in plots}

    def view(self) -> None:
        """Render the multiview selector, execute plot."""
        which = util.DictSelect(
            "View", self.plots,
            format_func=lambda x: self.plotnames[x],
            context=self.context)

        which.view()


class TermiteViewAltair(TermiteView):
    """View specifically for altair plots."""

    pass


@dataclass
class CountPlotAltair(TermiteViewAltair):
    """Barplot showing object counts."""

    name: str = "Count Plot (Altair)"

    def view(self) -> None:
        """Render the view."""
        col = self.get_catcol()

        colval = self.data[col].unique()
        sortorder = try_num_sort(colval)
        c = alt.Chart(self.data)\
               .mark_bar()\
               .encode(
                   alt.X(col, type='nominal', sort=sortorder),
                   alt.Y(col, type='quantitative', aggregate='count'),
        )

        self.ctxPlot.altair_chart(c, use_container_width=True)


@dataclass
class EcdfPlotAltair(TermiteViewAltair):
    """Altair based ecdf plot."""

    name: str = "ECDF Plot (Altair)"

    def view(self) -> None:
        """Render the plot."""
        col = self.get_numcol()

        c = alt.Chart(self.data).transform_window(
            ecdf="cume_dist()",
            sort=[{"field": col}],
        ).mark_line(
            interpolate="step-after"
        ).encode(
            alt.X(col, type="quantitative"),
            alt.Y("ecdf", type="quantitative"),
        )

        self.ctxPlot.altair_chart(c, use_container_width=True)


@dataclass
class ScatterPlotAltair(TermiteViewAltair):
    """Scatter plot (Altair).


    Not yet - too slow?
    """

    name: str = "Scatter Plot (Altair)"

    def view(self) -> None:
        """Render the plot."""
        cx = self.get_numcol(colname='x')
        cy = self.get_numcol(colname='y')

        c = alt.Chart(self.data)\
               .mark_point()\
               .encode(
                   x=cx,
                   y=cy,
               )

        self.ctxPlot.altair_chart(c, use_container_width=True)


@dataclass
class BoxPlotAltair(TermiteViewAltair):
    """Altair based boxplot."""

    name: str = "BoxPlot (Altair)"

    def view(self) -> None:
        """Render the plot."""
        numcol = self.get_numcol()
        catcol = self.get_catcol()

        colval = self.data[catcol].unique()
        sortorder = try_num_sort(colval)

        c = alt.Chart(self.data)\
               .mark_boxplot(extent="min-max")\
               .encode(
                   alt.X(numcol, type="quantitative").scale(zero=False),
                   alt.Y(catcol, type="nominal", sort=sortorder),
        )

        self.ctxPlot.altair_chart(c, use_container_width=True)


@dataclass
class ViolinPlotAltair(TermiteViewAltair):
    """Show an altair based violin plot."""

    name: str = "Violin Plot (Altair)"

    def view(self) -> None:
        """Render the plot."""
        numcol = self.get_numcol()
        catcol = self.get_catcol()

        colval = self.data[catcol].unique()
        sortorder = try_num_sort(colval)

        c = alt.Chart(self.data, height=50)\
               .transform_density(
                   numcol,
                   as_=[numcol, 'density'],
                   groupby=[catcol])\
               .mark_area(orient='vertical').encode(
                   alt.Y('density', type='quantitative')
            .stack('center')
            .impute(None)
            .title(None)
            .axis(labels=False, values=[0], grid=True, ticks=True),
                   alt.X(numcol, type='quantitative'),
                   alt.Row(catcol, type='nominal')
            .header(labelAlign='left',
                               labelAngle=0, )
        )

        self.ctxPlot.altair_chart(c, use_container_width=True)


@dataclass
class TopGenes(TermiteView):
    """Render barplot specifically for top expressing genes."""

    name: str = "Barplot of top expressed genes"

    def view(self) -> None:
        """Render the plot."""
        fig = px.bar(self.data, y='sumval', x='gene',
                     title="Most highly expressed genes")
        self.ctxPlot.plotly_chart(fig, use_container_width=True)


@dataclass
class Table(TermiteView):
    """Show a table."""

    name: str = "Table"

    def view(self) -> None:
        """Render the plot."""
        self.ctxPlot.dataframe(self.data)


@dataclass
class BarPlot(TermiteView):
    """Show a Barplot."""

    name: str = "Barplot"

    def view(self) -> None:
        """Render the plot."""
        fig = px.bar(self.data,
                     x=self.data.columns[0],
                     y=self.data.columns[1],
                     title=self.title,
                     text_auto='.3s')
        fig.update_xaxes(type='category')
        self.ctxPlot.plotly_chart(fig, use_container_width=True)


@dataclass
class DensityHeatmap(TermiteView):
    """Show a 2d heatmap of densities."""

    name: str = "Density Heatmap"

    def view(self) -> None:
        """Render the plot."""
        fig = px.density_heatmap(
            self.data, x=self.data.columns[0], y=self.data.columns[1],
            text_auto=True, nbinsx=20, nbinsy=20)
        self.ctxPlot.plotly_chart(fig, use_container_width=True)


@dataclass
class DensityHeatmap2(TermiteView):
    """Another density heatmap."""

    name: str = "Density Heatmap Alt"

    def view(self) -> None:
        """Render the plot."""
        c1, c2 = self.data.columns[:2]

        col1, col2 = self.ctxParam.columns(2)
        nobins = col1.slider("No Bins", min_value=10,
                             max_value=100, value=25, step=5)
        boxc = col2.container()
        logx = boxc.checkbox(f"Log X / {c2}")
        logy = boxc.checkbox(f"Log Y / {c1}")
        # logz = boxc.checkbox("Log Z")

        if logx:
            self.data[c2] = np.log1p(self.data[c2])
        if logy:
            self.data[c1] = np.log1p(self.data[c1])
        s1 = pd.cut(self.data.iloc[:, 0], bins=nobins,
                    retbins=False, duplicates='drop')
        s2 = pd.cut(self.data.iloc[:, 1], bins=nobins,
                    retbins=False, duplicates='drop')

        bd = pd.DataFrame(
            {c1: s1.apply(lambda x: x.left),
             c2: s2.apply(lambda x: x.left),
             }).pivot_table(index=c1, columns=c2, aggfunc=len)

        fig = px.imshow(bd, aspect='auto', height=800)

        self.ctxPlot.plotly_chart(fig, use_container_width=True)


@dataclass
class HeatMap(TermiteView):
    """Heatmap I guess."""

    name: str = "Heatmap"

    def view(self) -> None:
        """Render the plot."""
        data = self.data
        col1, col2 = self.ctxParam.columns(2)
        transpose = col2.checkbox("Transpose", value=False)

        norm1 = f"{data.index.name} (to %)"
        norm2 = f"{data.columns.name} (to %)"
        norm = str(col1.radio(
            "Normalize",
            ["None", norm1, norm2]
        ))

        if norm.startswith("None"):
            pass
        elif norm == norm1:
            data = 100 * data.divide(data.sum(1), axis=0)
        elif norm == norm2:
            data = 100 * data.divide(data.sum(0), axis=1)

        if transpose:
            data = data.T

        fig = px.imshow(data, text_auto=True, title=self.title)
        fig.update_xaxes(type='category')
        fig.update_yaxes(type='category')

        self.ctxPlot.plotly_chart(fig, use_container_width=True)


@dataclass
class Ecdf(TermiteView):
    """Another ecdf plot."""

    name: str = "Ecdf plot"
    usecol: int = 0

    def view(self) -> None:
        """Render the plot."""
        rugplot = False
        showrug = 'rug' if rugplot else None

        data = self.data.columns[self.usecol]
        fig = px.ecdf(self.data, title=self.title,
                      marginal=showrug, ecdfnorm='percent',
                      x=data)
        fig.update_traces(hoverinfo='none')
        self.ctxPlot.plotly_chart(fig, use_container_width=True)


@dataclass
class SnsEcdf(TermiteView):
    """Show a seaborn ecdf plot."""

    name: str = "Ecdf plot"
    usecol: int = 0

    def view(self) -> None:
        """Render the plot."""
        ax = sns.ecdfplot(self.data)
        st.pyplot(ax.figure)


@dataclass
class Histogram(TermiteView):
    """Show histogram."""

    name: str = "Histogram"

    def view(self) -> None:
        """Render the plot."""
        nozero = int((self.data == 0).sum())
        datalen = len(self.data)
        skip0 = self.ctxParam.checkbox(
            f"Hide zero ({nozero} / {datalen} = {100*nozero/datalen:.1f}%)",
            value=False, key='hide0')
        log_y = self.ctxParam.checkbox(
            "Log y axis", value=False, key='logy')

        data = self.data.copy()
        if skip0:
            data = data[data.iloc[:, 0] != 0]

        fig = px.histogram(data, log_y=log_y, x=self.data.columns[0])
        self.ctxPlot.plotly_chart(fig, use_container_width=True)


class HeatmapAltair(TermiteViewAltair):
    """Altair based heatmap."""

    name: str = "Heatmap (altair)"

    def view(self) -> None:
        """Render the plot."""
        pass


@dataclass
class Scatter(TermiteView):
    """Show a scatter plot."""

    name: str = "Scatter"
    pointsize: int = 2
    zloggable: bool = False
    xloggable: bool = True
    yloggable: bool = True

    def view(self) -> None:
        """Render the plot."""
        cx = self.get_numcol(colname='x')
        cy = self.get_numcol(colname='y', exclude=[cx])

        col_options = self.numcol + self.catcol

        def _color_sort(c: Any) -> str:
            """Sort keys to force cx, cy to the end."""
            if c in [cx, cy]:
                return 'zzzzzz' + str(c)
            return str(c)
        # move cx, cy to the end
        col_options = sorted(col_options, key=_color_sort)

        col = util.selectbox_mem(
            context=self.ctxParam1,
            label='Color on', key='color',
            options=col_options)

        point_size = self.ctxParam1.slider(
            "Point size", min_value=1, max_value=20, value=self.pointsize,
            step=2, key='pointsize')
        opacity = self.ctxParam1.slider(
            "Opacity", min_value=0.1, max_value=1.0, value=0.8,
            step=0.1, key='opacity')
        ccmap = util.selectbox_mem(
            context=self.ctxParam1,
            label="Cont. Colorscale",
            options=['viridis_r', 'Plotly3_r', 'Cividis_r', 'Agsunset_r'],
            key='ccmap')

        c1, c2, c3 = self.ctxParam1.columns(3)

        if self.xloggable \
           and util.checkbox_mem(
               context=c1, label=f'Log X ({cx})', key='logx'):
            self.data[cx] = np.log1p(self.data[cx])

        if self.yloggable and \
           util.checkbox_mem(
               context=c2, label=f'Log Y ({cy})', key='logy'):
            self.data[cy] = np.log1p(self.data[cy])

        if self.zloggable:
            if util.checkbox_mem(
                    context=c3, label=f'Log Z ({col})', key='col'):
                self.data[col] = np.log1p(self.data[col])

        fig = px.scatter(self.data,
                         x=cx, y=cy, color=col,
                         color_continuous_scale=ccmap,
                         color_discrete_sequence=px.colors.qualitative.Dark24,
                         height=800)

        fig.update_traces(marker=dict(size=point_size, opacity=opacity))

        self.ctxPlot.plotly_chart(fig, use_container_width=True)
