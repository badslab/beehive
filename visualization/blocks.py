
from functools import partial
from textwrap import dedent
from typing import Tuple, Optional, Generator
from urllib.parse import quote


from numpy import False_
import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from termite.vis import util
from termite.vis.util import persist
from termite.vis import data as vdata


def get_dataset(with_search_string: bool = False
                ) -> pd.Series:
    if with_search_string:
        search_string = st.sidebar.text_input("Search string")
        experiment = util.SqlSelect(
            "Experiment", key='exp',)
    else:
        col1, col2 = st.sidebar.columns([7,3])
        col1.caption('**Select Experiment/Dataset**')
        if col2.toggle(label=":mag:"):
            col1, col2 = st.sidebar.columns([3,7])
            col1.caption("Filter:")
            search_string = col2.text_input(
                "Search string", label_visibility='collapsed')
            sql= f"""SELECT DISTINCT experiment, author, title, experiment
                       FROM dataset_md
                      WHERE ( title ILIKE '%{search_string}%'
                         OR author ILIKE '%{search_string}%'
                         OR experiment ILIKE '%{search_string}%' )
                       LIMIT 100 """
        else:
            sql= """SELECT DISTINCT experiment, author, title, experiment
                      FROM dataset_md"""

        col1, col2 = st.sidebar.columns([3,7])
        col1.caption('Experiment:')
        experiment = util.SqlSelect(
            "Experiment", key='exp', context=col2,
            label_visibility='collapsed', sql=sql)

    col1, col2 = st.sidebar.columns([3,7])
    col1.caption("Dataset:")
    dataset = util.SqlSelect(
        "Dataset", key='dataset',
        context=col2, label_visibility='collapsed',
        oformat="{layertype} ({layername})",
        sql=f"""SELECT dataset, layername, layertype,
                       starts_with(layertype, 'raw') as israw
                  FROM dataset_md
                 WHERE experiment='{experiment}'
              ORDER BY israw, layername
                 LIMIT 100 """)
    dataset_rec = util.SqlOneRecord(
        sql=f"""SELECT *
                  FROM dataset_md
                 WHERE dataset='{dataset}' """)
    title(dataset_rec)
    return dataset_rec


def title(dataset_rec: pd.Series,
          context: DeltaGenerator = st):

    ds = dataset_rec.loc['dataset']
    exp =  dataset_rec.loc['experiment']
    ds_quote = quote(ds)
    exp_quote = quote(exp)
    short_title = util.short_title(dataset_rec['title'])
    year = dataset_rec.loc['year']
    ds_url = f"/?app=dataset&dataset={ds_quote}&exp={exp_quote}"

    context.markdown(dedent(f"""
         <a href="{ds_url}" target="_self">
           <b>{short_title}</b>
           <i>{dataset_rec['author']}</i>
           ({year}, {dataset_rec['dataset']})</a>"""),
        unsafe_allow_html=True)


def get_gcn(name: str,
            ds_id: str,
            exp_id: str,
            key: str = "",
            context: DeltaGenerator = st.sidebar,
            allow_dimred: bool = False,
            allow_nothing: bool = False,
            ) -> Generator[Tuple[str, str, pd.Series], None, None]:
    """Return a gene, numer. categ. or dimred column(s).

    Be careful with dimred's - these yield two columns!
    """
    ctype_name = dict(
        none='None',
        cat='Categorical',
        num='Numerical',
        gen='Gene',
        dimred='Dim. Red.')

    if not allow_nothing:
        del ctype_name['none']
    if not allow_dimred:
        del ctype_name['dimred']

    col1, col2, col3 = context.columns([3, 4, 3])
    col1.caption(f"**{name}**")
    ctype = util.selectbox_mem(
        label='what', options=ctype_name.keys(),
        format_func=lambda x: ctype_name[x],
        key=f"{key}_t", context=col2,
        label_visibility='collapsed')
    if ctype != 'dimred':
        force_filter: bool \
            = col3.toggle(label=":tornado:", key=f"{key}_f")
    else:
        force_filter = False

    if ctype == "none":
        yield ctype, "", pd.Series()
        return

    if key:
        key += '_'

    gdargs = dict(
        force_filter=force_filter,
        hide_caption=True)

    if ctype == 'cat':
        yield ctype, *get_categorical(
            exp_id, key=f"{key}_cat", **gdargs)
    elif ctype == 'num':
        yield ctype, *get_numerical(
            exp_id, key=f"{key}_num", **gdargs)
    elif ctype == 'gen':
        yield ctype, *get_expr(
            ds_id, key=f"{key}_gene", **gdargs)
    elif ctype == 'dimred':
        _, col1, col2 = st.sidebar.columns([5, 25, 70])
        col1.caption('Name:')
        dimred = get_dimred_name(
            exp_id=exp_id,
            context=col2,
            label_visibility='collapsed')

        dim1 = f"{dimred}/00"
        dim2 = f"{dimred}/01"

        yield 'dimred', dim1, get_obs_num(exp_id=exp_id, colname=dim1)
        yield 'dimred', dim2, get_obs_num(exp_id=exp_id, colname=dim2)


def colquery(nfunc, dfunc, key: str,
             force_filter: Optional[bool] = None):
    """
    Perform a generic query on a column.

    Relies on two curried functions:
     * nfunc - queries for the name of the variable to retrieve
     * dfunc - gets the actual data.
    """

    if force_filter is None:
        col1, col2 = st.sidebar.columns([85, 15])
        use_filter = col2.button(":mag:", key=f"use_filter_{key}")
        data_name = nfunc(name=key, key=key, context=col1)
    else:
        _, col1, col2 = st.sidebar.columns([0.5, 2.5, 7])
        col1.caption('Name:')
        data_name = nfunc(name=key, key=key, context=col2)

    # get the data
    data = dfunc(name=data_name)

    # do we need to filter?
    fkey = f"filter_{key}"
    qp = st.experimental_get_query_params()

    if (force_filter is None and (use_filter or fkey in qp)) or force_filter:
        # Get a filter to query on - see the pandas query function
        if force_filter is None:
            col1, col2 = st.sidebar.columns([30, 70])
        else:
            _, col1, col2 = st.sidebar.columns([5, 25, 70])
        col1.caption('Filter:')
        filter_str = persist(
            col2.text_input,
            label=fkey, key=fkey,
            label_visibility='collapsed')
        filter_str = filter_str.strip()

        if filter_str:
            # we actually have a filter!
            cc = pd.DataFrame(dict(x=data))
            try:
                cc = cc.query(filter_str)
            except Exception as e:
                st.error(dedent(
                    f"""
                    Invalid filter: {filter_str}\n\n
                    {type(e).__name__} : {e}
                    Suggestions:
                        x == "something"
                        x > 100 and x < 200
                    """))
                st.stop()
            data = cc['x']

    # be sure  that the returned series has the correct name
    data.name = data_name
    return data_name, data


def get_expr(ds_id: str,
             name: str = "Gene",
             key: str ="gene",
             force_filter: bool = False,
             hide_caption:bool = False) \
             -> Tuple[str, pd.Series]:
    """Get expression data of a gene."""
    if not hide_caption:
        st.sidebar.caption(name)

    sql = f"""SELECT DISTINCT gene
                FROM help_gene
               WHERE dataset_id='{ds_id}'
                ORDER BY gene"""

    return colquery(
        nfunc = partial(util.SqlSelect, sql=sql, label_visibility="collapsed"),
        dfunc = partial(vdata.get_expr, ds_id=ds_id),
        force_filter=force_filter,
        key=key)


def get_categorical(
        exp_id: str,
        name: str = "Categorical",
        key: str = "cat",
        hide_caption: bool = False,
        force_filter: bool = False,
        exclude: list = [] ) -> Tuple[str, pd.Series]:

    # prep sql exclude statement to not show those
    # in the dropdown
    if len(exclude) > 0:
        join = "','".join(exclude)
        exclude_stmt = f"AND name NOT IN ('{join}')"
    else:
        exclude_stmt = ""

    sql = f"""SELECT name, name
                  FROM help_obs
                 WHERE exp_id='{exp_id}'
                   AND dtype='categorical'
                       {exclude_stmt} """

    if not hide_caption:
        st.sidebar.caption(name)

    return colquery(
        nfunc = partial(util.SqlSelect, sql=sql, label_visibility="collapsed"),
        dfunc = partial(vdata.get_obs_cat, exp_id=exp_id),
        force_filter=force_filter,
        key=key)


def get_numerical(
        exp_id: str,
        name: str = "Numerical",
        key: str = "num",
        hide_caption: bool = False,
        force_filter: bool = False,
        exclude: list[str] = []) -> Tuple[str, pd.Series]:


    # prep sql exclude statement to not show those
    # in the dropdown
    if len(exclude) > 0:
        join = "','".join(exclude)
        exclude_stmt = f"AND name NOT IN ('{join}')"
    else:
        exclude_stmt = ""

    sql=f"""SELECT name, name
              FROM help_obs
             WHERE exp_id='{exp_id}'
               AND dtype IN ('int', 'float')
                   {exclude_stmt} """

    if not hide_caption:
        st.sidebar.caption(name)

    return colquery(
        nfunc = partial(util.SqlSelect, sql=sql, label_visibility="collapsed"),
        dfunc = partial(vdata.get_obs_num, exp_id=exp_id),
        force_filter=force_filter,
        key = key)


def get_expr2(ds_id: str,
              gene1: str,
              gene2: str) -> pd.DataFrame:
    """Return expression of two genes."""
    return util.execute_sql(
        sql=f"""SELECT a.value as g1, b.value as g2
                  FROM expr as a
                  JOIN expr as b on (a.obs=b.obs)
                  WHERE a.dataset_id={ds_id}
                    AND b.dataset_id={ds_id}
                    AND a.gene = '{gene1}'
                    AND b.gene='{gene2}'
             """)


def get_obs_num(exp_id: str,
                colname: str) -> pd.Series:
    """Return an numerical obs column."""
    rv = util.execute_sql(
        sql=f"""SELECT obs, value
                  FROM obs_num
                 WHERE exp_id={exp_id}
                   AND name='{colname}' """)
    rvcol = rv.set_index('obs')['value']
    rvcol.name = colname
    return rvcol


def get_dimred_name(
        exp_id: str,
        context: DeltaGenerator = st.sidebar,
        label_visibility: str = 'visible') -> str:
    """Return one dim.red name."""
    rv = util.SqlSelect(
        "Dimred", context=context,
        label_visibility=label_visibility,
        sql=f"""
            SELECT DISTINCT dimred_name
              FROM help_obs
             WHERE exp_id={exp_id}
               AND dtype='dimred' """)
    assert isinstance(rv, str)
    return rv
