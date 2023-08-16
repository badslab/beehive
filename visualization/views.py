
from typing import Optional
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from termite.vis import util


@dataclass
class TermiteView:
    data: pd.DataFrame
    name: str = "Termite plot"
    title: Optional[str] = None
    subject: Optional[str] = None
    
    ctxPlot: DeltaGenerator = field(default_factory=st.empty)
    ctxData: DeltaGenerator = field(default_factory=st.container)
    ctxParam: DeltaGenerator = field(default_factory=st.container)

    def __post_init__(self):
        if self.ctxPlot is None:
            self.ctxPlot = st.empty()
        if self.ctxData is None:
            self.ctxData = st.container()
        if self.ctxParam is None:
            self.ctxParam = st.container()
        if self.title is None:
            title = self.name
            if self.subject is not None:
                title += f" of {self.subject}"
            self.title = title
            
class Multi():
    def __init__(self, *plots):
        self.plots = {a.name: a for a in plots}
        
    def view(self):
        which = util.DictSelect("View", self.plots)
        which.view()


@dataclass
class TopGenes(TermiteView):
    name: str = "Barplot of top expressed genes"
    
    def view(self):
        fig = px.bar(self.data, y='sumval', x='gene', 
                     title="Most highly expressed genes")
        self.ctxPlot.plotly_chart(fig, use_container_width=True)

        
@dataclass   
class BarPlot(TermiteView):
    name: str = "Barplot"
    def view(self):
        fig = px.bar(self.data,
                     x=self.data.columns[0],
                     y=self.data.columns[1],
                     title=self.title,
                     text_auto='.3s')
        fig.update_xaxes(type='category')
        self.ctxPlot.plotly_chart(fig, use_container_width=True)

@dataclass
class DensityHeatmap(TermiteView):
    name: str = "Density Heatmap"

    def view(self):
        
        fig = px.density_heatmap(
            self.data, x=self.data.columns[0], y=self.data.columns[1],
            text_auto=True, nbinsx=20, nbinsy=20)
        self.ctxPlot.plotly_chart(fig, use_container_width=True)

        
@dataclass
class DensityHeatmap2(TermiteView):
    name: str = "Density Heatmap Alt"

    def view(self):
        
        c1, c2 = self.data.columns[:2]

        col1, col2 = self.ctxParam.columns(2)
        nobins = col1.slider("No Bins", min_value=10, max_value=100, value=25, step=5)
        boxc = col2.container()
        logx = boxc.checkbox(f"Log X / {c2}")            
        logy = boxc.checkbox(f"Log Y / {c1}")
        logz = boxc.checkbox("Log Z")
            
        if logx:
            self.data[c2] = np.log1p(self.data[c2])
        if logy:
            self.data[c1] = np.log1p(self.data[c1])
        s1 = pd.cut(self.data.iloc[:,0], bins = nobins, retbins=False, duplicates='drop')
        s2 = pd.cut(self.data.iloc[:,1], bins = nobins, retbins=False, duplicates='drop')

        bd = pd.DataFrame(
            {c1:s1.apply(lambda x: x.left),
             c2:s2.apply(lambda x: x.left),
             }).pivot_table(index=c1, columns=c2, aggfunc=len)

        
        fig = px.imshow(bd, aspect='auto', height=800)

        self.ctxPlot.plotly_chart(fig, use_container_width=True)
        
@dataclass
class HeatMap(TermiteView):
    name: str = "Heatmap"

    def view(self):
        data = self.data
        col1, col2 = self.ctxParam.columns(2)
        norm = str(col1.radio(
            "Normalize",
            [f"None",
             f"Rows ({data.index.name}) to %",
             f"Columns ({data.columns.name}) to %"]))
        
        transpose = col2.checkbox("Transpose", value=False)

        if norm.startswith("None"):
            pass
        elif norm.startswith("Rows"):
            data = 100 * data.divide(data.sum(1), axis=0) 
        elif norm.startswith("Columns"):
            data = 100 * data.divide(data.sum(0), axis=1)

        if transpose:
            data = data.T
            
        fig = px.imshow(data, text_auto=True, title=self.title)
        fig.update_xaxes(type='category')
        fig.update_yaxes(type='category')

        self.ctxPlot.plotly_chart(fig, use_container_width=True)

@dataclass
class Ecdf(TermiteView):
    name: str = "Ecdf plot"
    usecol: int = 0
    
    def view(self):
        rugplot = self.ctxParam.checkbox(
            "Show Rugplot", value=False, key='rugplot')
        showrug = 'rug' if rugplot else None
        if isinstance(self.data, pd.Series):
            data = self.data
        else:
            data = self.data.columns[self.usecol]
        fig = px.ecdf(self.data, title=self.title,
                      marginal=showrug, ecdfnorm='percent',
                      x=data)
        self.ctxPlot.plotly_chart(fig, use_container_width=True)

        
@dataclass
class Histogram(TermiteView):
    name: str = "Histogram"
    def view(self):
        skip0 = self.ctxParam.checkbox(
            "Hide zero", value=False, key='hide0')
        log_y = self.ctxParam.checkbox(
            "Log y axis", value=False, key='logy')
        
        data = self.data.copy()
        if skip0:
            data = data[data.iloc[:,0] != 0]

        fig = px.histogram(data, log_y=log_y, x=self.data.columns[0])
        self.ctxPlot.plotly_chart(fig, use_container_width=True)
        

@dataclass
class Table(TermiteView):
    def view(self):
        self.ctxPlot.dataframe(
            self.data, use_container_width=True)
    

def scatter(data,           
            ctxPlot: DeltaGenerator,
            ctxParam: DeltaGenerator,
            labels=None):

    c1, c2 = ctxParam.columns(2)
    
    point_size = c1.slider("Point size", min_value=1, max_value=20,
                                 value=4, step=2, key='pointsize')
    opacity = c2.slider("Opacity", min_value=0.1, max_value=1.0,
                                 value=0.8, step=0.1, key='opacity')
    
    fig = px.scatter(data,
                     x=data.columns[0],
                     y=data.columns[1],
                     color=data.columns[2],
                     labels=labels,
                     height=800)

    fig.update_traces(marker=dict(size=point_size, opacity=opacity))
    
    ctxPlot.plotly_chart(fig, use_container_width=True)


