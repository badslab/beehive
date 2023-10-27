"""Utilities for Streamlit Visualization."""

import logging
from typing import Any, Dict, Optional, Callable

from duckdb import DuckDBPyConnection
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from termite.db import get_conn

lg = logging.getLogger(__name__)
lg.setLevel(logging.DEBUG)

COMMONWORDS = """the at there some my of be use her than and this an
 would first a have each make water to from which like been in or she
 him call is one do into who you had how time oil that by their has
 its it word if look now he but will two find was not up more long for
 what other write down on all about go day are were out see did as we
 many number get with when then no come his your them way made they
 can these could may I said so people part """.split()


def short_title(title: str) -> str:
    """Create a short title (from a longer one)."""
    short_title = " ".join(title[:80].split()[:-1])
    while True:
        a, b = short_title.rsplit(None, 1)
        if b in COMMONWORDS:
            short_title = a
        else:
            break
    if short_title != title:
        short_title += '...'
    return short_title


def item_title(item_name: str, item: Dict[str, Any]) -> str:
    """Return title, otherwise make one."""
    if 'title' in item:
        return str(item['title'])
    else:
        return item_name.replace('_', ' ').capitalize()


# streamlit components that get/copy state from/to the url bar
# making bookmarkable views..
def selectbox_mem(
        context: DeltaGenerator,
        label: str,
        options: list[str],
        label_visibility: str = "visible",
        default: Optional[str] = None,
        format_func: Optional[Callable] = None,
        index: int = 0,
        key: Optional[str] = None) -> str:
    """URL persistent selectbox."""
    if key is None:
        key = label.lower().replace(' ', '_')

    options = list(options)
    qp = st.experimental_get_query_params()

    def update_query_param() -> None:
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    # find out what index to point to
    idx = 0
    if key in qp:
        # do we have the key in the URL?
        # url takes precedence
        qval = qp[key][0]
        if qval in options:
            idx = options.index(qval)
    elif default is not None and default in options:
        # is there a default specified?
        idx = options.index(default)
    else:
        idx = index

    sbargs = {}
    if format_func is not None:
        sbargs['format_func'] = format_func

    return context.selectbox(label, options, key=key, index=idx,
                             label_visibility=label_visibility,
                             on_change=update_query_param, **sbargs)


def persist(
        widget, label,
        default=None,
        **kwargs):

    kwargs['key'] = key = kwargs.get(
        'key', label.lower().replace(' ', ''))

    qp = st.experimental_get_query_params()

    def update_query_param() -> None:
        """Update changed value in query parameters."""
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    kwargs['on_change'] = update_query_param

    if widget.__name__ == 'text_input':
        convert_from_qp = str
        def convert_default(x): return "" if x is None else str(x)

    if key in qp:
        kwargs['value'] = convert_from_qp(qp[key][0])
    else:
        kwargs['value'] = convert_default(default)

    return widget(label=label,
                  **kwargs)


def textbox_mem(
        context: DeltaGenerator,
        label: str,
        default: str = "",
        key: Optional[str] = None) -> str:

    if key is None:
        key = label.lower().replace(' ', '_')

    qp = st.experimental_get_query_params()

    def update_query_param() -> None:
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    # set default value via session state - otherwise we get
    # errors when we update this through other methods
    if key in qp:
        dvalue = qp[key][0]
    else:
        dvalue = default

    return context.text_input(label, key=key, value=dvalue,
                              on_change=update_query_param)


def checkbox_mem(
        context: DeltaGenerator,
        label: str,
        default: bool = False,
        key: Optional[str] = None) -> bool:
    """checkbox that persists in the URL."""
    if key is None:
        key = label.lower().replace(' ', '')

    qp = st.experimental_get_query_params()

    def update_query_param() -> None:
        """Set changed value in the query parameters."""
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    # set default value via session state - otherwise we get
    # errors when we update this through other methods
    dvalue = default
    if key in qp:
        dvalue = bool(qp[key][0])

    return context.checkbox(label, key=key, value=dvalue,
                            on_change=update_query_param)


def execute_sql(
        sql: str,
        db: Optional[DuckDBPyConnection] = None) -> pd.DataFrame:
    """Execute SQL, return a pandas dataframe."""

    if db is None:
        db = get_conn()

    lg.debug(f"executing sql:\n{sql}\n--------------")
    df = db.execute(sql).df()
    lg.debug(f"result: {len(df)} records")
    return df


def SqlOneRecord(sql) -> pd.Series:
    """Return just one record - error if none or more."""
    rv = execute_sql(sql)
    if len(rv) == 0:
        st.error("No records found for {sql}")
        st.stop()
    elif len(rv) > 1:
        st.error("More than 1 records found for {sql}")
        st.write(rv)
        st.stop()
    else:
        return rv.iloc[0]


def DictSelect(name, data, format_func=None,
               key=None, context: DeltaGenerator = st.sidebar,
               label_visibility: str = "visible", ) -> Any:
    """Show widget, based on dict."""
    if key is None:
        key = "_".join(name.lower().split())
    options = list(data.keys())
    key = selectbox_mem(
        context, name, options, key=key,
        format_func=format_func,
        label_visibility=label_visibility)
    return data[key]


def SqlSelect(name, sql, key=None,
              context: DeltaGenerator = st.sidebar,
              label_visibility: str = "visible",
              oformat: Optional[str] = None) -> Any:
    """Show selectbox based on sql query."""
    df = execute_sql(sql)
    if len(df) == 0:
        st.warning("No hits found")
        st.stop()

    sbargs = dict()

    if len(df.columns) == 1:
        options = list(df.iloc[:, 0])
    else:
        if oformat is None:
            options_dict: Dict[str, str] = {
                r.iloc[0]: ' / '.join(map(str, r.iloc[1:]))
                for (_, r) in df.iterrows()}
        else:
            options_dict: Dict[str, str] = {
                r.iloc[0]: oformat.format(**r)
                for (_, r) in df.iterrows()}
        options = list(options_dict.keys())
        sbargs['format_func'] = lambda x: options_dict[x]

    return selectbox_mem(
        context, name, options,
        label_visibility=label_visibility,
        key=key, **sbargs)
