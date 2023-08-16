
import logging
from typing import Any, Dict, Optional, Callable

from duckdb import DuckDBPyConnection
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

lg = logging.getLogger(__name__)
lg.setLevel(logging.DEBUG)


from termite.db import get_conn

COMMONWORDS = """the at there some my of be use her than and this an
 would first a have each make water to from which like been in or she
 him call is one do into who you had how time oil that by their has
 its it word if look now he but will two find was not up more long for
 what other write down on all about go day are were out see did as we
 many number get with when then no come his your them way made they
 can these could may I said so people part """.split()


def short_title(title):
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

def selectbox_mem(
        context: DeltaGenerator,
        label: str,
        options: list[str],
        default: Optional[str] = None,
        format_func: Optional[Callable] = None,
        index: int =0,
        key: Optional[str] = None) -> Any:

    if key is None:
        key = label.lower().replace(' ', '_')
        
    options = list(options)
    qp = st.experimental_get_query_params()

    def update_query_param():
        qp[key] = st.session_state[key]
        st.experimental_set_query_params(**qp)

    if default is not None and default in options:
        idx = options.index(default)
    else:
        idx = index
        
    if key in qp:
        qval = qp[key][0]
        if qval in options:
            idx = options.index(qval)

    sbargs = {}
    if format_func is not None:
        sbargs['format_func'] = format_func
        
    return context.selectbox(label, options, key=key, index=idx,
                             on_change=update_query_param, **sbargs)


def item_title(item_name, item):
    if 'title' in item:
        return item['title']
    else:
        return item_name.replace('_', ' ').capitalize()

    
def textbox_mem(
        context: DeltaGenerator,
        label: str,
        default: str="",
        key: Optional[str]=None) -> str:

    if key is None:
        key = label.lower().replace(' ', '_')
        
    qp = st.experimental_get_query_params()

    def update_query_param():
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


def execute_sql(
        sql: str,
        db: Optional[DuckDBPyConnection] = None) -> pd.DataFrame:
    
    """Execute SQL, return a pandas dataframe"""

    if db is None:
        db = get_conn()
        
    lg.debug(f"executing sql:\n{sql}\n--------------")
    df = db.execute(sql).df()
    lg.debug(f"result: {len(df)} records")
    return df


def SqlOneRecord(sql):
    """Return just one record - error if none or more"""
    rv = execute_sql(sql)
    if len(rv) == 0:
        st.error("No records found for {sql}")
        st.stop()
    elif len(rv) > 1:
        st.error("More than 1 records found for {sql}")
        st.write(rv)
        st.stop()
    else: return rv.iloc[0]
    

def DictSelect(name, data, format_func=None, key=None):
    if key is None:
        key = "_".join(name.lower().split())
    options = list(data.keys())
    key = selectbox_mem(
        st.sidebar, name, options, key=key,
        format_func=format_func)
    return data[key]


def SqlSelect(name, sql, key=None,
              oformat:Optional[str] = None) -> Any:
        
    df = execute_sql(sql)
    if len(df) == 0:
        st.warning("No hits found")
        st.stop() 
    
    sbargs = dict()

    if len(df.columns) == 1:
        options = list(df.iloc[:,0])
    else:
        if oformat is None:
            options_dict: Dict[str, str] = {
                r.iloc[0]:' / '.join(map(str, r.iloc[1:]))
                for (_, r) in df.iterrows()}
        else:
            options_dict: Dict[str, str] = {
                r.iloc[0]:oformat.format(**r)
                for (_, r) in df.iterrows()}
        options = list(options_dict.keys())
        sbargs['format_func'] = lambda x: options_dict[x]

    return selectbox_mem(
        st.sidebar, name, options, key=key, **sbargs)

