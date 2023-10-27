"""Generic table parser & viewer."""

from collections import defaultdict
import logging
from typing import List

import pandas as pd
import streamlit as st

import views

import blocks

lg = logging.getLogger(__name__)


def gentab() -> None:
    """Streamlit generic table viewer function."""
    dataset_rec = blocks.get_dataset()
    ds_id = dataset_rec.loc['dataset_id']
    exp_id = dataset_rec.loc['experiment_id']

    viewbox = st.sidebar.empty()

    st.sidebar.markdown('**Data Columns**')
    columns = defaultdict(list)
    data: pd.DataFrame

    done_columns = False
    colno = 0
    while True:
        for ctype, cname, ccol in \
                blocks.get_gcn(
                    name=f"Column {colno}",
                    allow_dimred=True,
                    ds_id=ds_id,
                    exp_id=exp_id,
                    key=f"c{colno}",
                    allow_nothing=(colno > 1),
                ):

            colno += 1
            if ctype == "none":
                done_columns = True
                break

            columns[ctype].append(cname)
            if colno == 1:
                data = pd.DataFrame(ccol)
            else:
                data = data.merge(
                    ccol, how='inner',
                    left_index=True, right_index=True)

        if done_columns:
            break

    composition = {x: len(columns[x]) for x in 'cat num gen dimred'.split()}
    compstr = " ".join([f"{a}:{b}" for (a, b) in composition.items()])
    st.markdown(f"**Data points**: {data.shape[0]:,d} -- {compstr}")

    allviews: List[views.TermiteView] \
        = [views.Table(data=data)]

    nocat = composition['cat']
    nonum = composition['num'] + composition['gen'] + composition['dimred']

    if nocat > 0:
        allviews.append(
            views.CountPlotAltair(data=data)
        )

    if nonum > 0:
        allviews.append(
            views.EcdfPlotAltair(data=data)
        )

    if nonum > 1:
        allviews.append(
            views.Scatter(data=data)
        )

    if nocat > 0 and nonum > 0:
        allviews.extend([
            views.BoxPlotAltair(data=data),
            views.ViolinPlotAltair(data=data)
            ])

    views.Multi(
        *allviews,
        context=viewbox
    ).view()
