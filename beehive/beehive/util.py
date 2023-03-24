"""Helper functions for beehive."""

import hashlib
import logging
import sqlite3
import subprocess as sp
import time
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from bokeh.plotting import curdoc

import beehive

lg = logging.getLogger(__name__)
lg.setLevel(logging.DEBUG)


def get_geneset_folder() -> Path:
    """Return folder with the raw geneset database
    """

    geneset_folder = beehive.BASEDIR / 'geneset'
    if not geneset_folder.exists():
        geneset_folder.mkdir(parents=True)
    return geneset_folder


@lru_cache(16)
def get_geneset_db(
        dsid: str) -> sqlite3.Connection:
    """Return database with genesets
    """
    lg.info("opening a new database")
    geneset_db_folder = get_geneset_folder() / "gsdb"
    # geneset_db_folder = get_geneset_folder() / "db"

    if not geneset_db_folder.exists():
        geneset_db_folder.mkdir(parents=True)

    gsdb = geneset_db_folder / f'geneset_db_{dsid}.sqlite'
    lg.debug(f"Get geneset db {gsdb}")

    return sqlite3.connect(gsdb)


def dict_set(data, *args, overwrite=False):
    """Careful dict value setter
        - create parents if required
        - never overwrite a value if it is already there
          (unless overwrite is True)

    args should at least be of length 2 - all but the last are interpreted as
    keys, the last value of *args is the value to set

    so dict_set(yaml, key1, key2, val) sets yaml[key1][key2] = val

    """

    assert len(args) >= 2

    if len(args) > 2:
        # first organize partent keys
        parent, *rest = args
        if parent not in data:
            data[parent] = {}
        return dict_set(data[parent], *rest, overwrite=overwrite)

    key, val = args
    if overwrite or key not in data:
        data[key] = val
    else:
        # not setting - value exists
        pass


def timer(func):
    def wrapper(*arg, **kw):
        '''source: http://www.daniweb.com/code/snippet368.html'''
        t1 = time.time()
        res = func(*arg, **kw)
        t2 = time.time()
        lg.info(f'Function run {1000*(t2-t1):.2f}ms.')
        return res
    return wrapper


def find_prq(dsid, ptype, check_exists=True):
    """ Find parquet file based on dataset id
        and parquet type"""

    assert ptype in ['X', 'obs', 'var', 'gsea']
    name = f"{dsid}.{ptype}.prq"
    prq_dir = beehive.BASEDIR / 'prq'
    prq_file = prq_dir / name

    if prq_file.exists():
        # lg.debug(f"Found PRQ file in prq/ {prq_file}")
        return prq_file

    # find alternative location:
    prq_file_alt = beehive.BASEDIR / 'data' / 'h5ad' / name
    if (not prq_file_alt.exists()) and check_exists:
        raise FileNotFoundError(f'Cannot find prq file {name}')

    lg.warning(f"Found Parquet file in data/h5ad {prq_file}")
    return prq_file


def get_datadir(name):
    """Return a bokeh view's data folder."""
    datadir = beehive.DATADIR / name
    lg.debug(f"Beehive data dir {datadir}")
    return datadir


def getarg(args, name, default=None, dtype=str):
    """Bokeh helper to get arguments from the URL."""

    if name in args:
        val = args[name][0].decode()
        # if name == "subset_categories":
        #     return
        return dtype(val)
    else:
        return default


#
# Bokeh helper functions
#
def list2options(optionlist, add_none=True):
    return [('-', '-')] + [(x, x) for x in optionlist]


def create_widget(name: str,
                  widget,
                  default=None,
                  title: Optional[str] = None,
                  curdoc=None,
                  update_url: bool = True,
                  value_type=None,
                  **kwargs):
    """Widget helper.

    Helper function that:
      - Creates a widget
      - Tries to find default value as a URL parameter, otherwise use `default`
      - Binds a javascript callback that update the URL bar upon change.

    Note - the js callback relies on a function called insertUrlParam
           being available.
    """
    import random

    from bokeh.models import RadioGroup  # type: ignore[attr-defined]
    from bokeh.models.callbacks import CustomJS
    from bokeh.models.widgets.inputs import AutocompleteInput

    assert curdoc is not None

    param_name = name  # to get & retrieve from the URL

    if title is None:
        title = name.capitalize()

    if widget == AutocompleteInput:
        # unique name - prevents browser based autocomplete
        # which messes up different datasets
        name = f"{name}.{random.randint(100000, 999999)}"

    js_onchange_value = CustomJS(
        code=f'insertUrlParam("{param_name}", cb_obj.value);')
    args = curdoc.session_context.request.arguments

    js_onchange_active = CustomJS(
        code=f'insertUrlParam("{param_name}", cb_obj.active);')

    if widget == RadioGroup:
        new_widget = widget(name=name, **kwargs)
        new_widget.active = getarg(args, param_name, default, dtype=value_type)

        if update_url:
            new_widget.js_on_change("active", js_onchange_active)
        return new_widget

    else:
        new_widget = widget(name=name, title=title, **kwargs)

    if (value_type == list) and not \
            (type(getarg(args, param_name, default)) == list):
        new_widget.value = getarg(args, param_name, default).split(",")
        if new_widget.value == [""]:
            new_widget.value = []
    elif (value_type == int) or (value_type == float):
        new_widget.value = getarg(args, param_name, default, dtype=value_type)
    else:
        new_widget.value = getarg(args, param_name, default)
    if update_url:
        new_widget.js_on_change("value", js_onchange_value)
    return new_widget

def get_dataset2():
    args = curdoc().session_context.request.arguments
    print(args)


def getcolor(x, palette, vmax, vmin=None):
    if vmin is None:
        vmin = -vmax

    v = int(255 * (x + vmax) / (2 * vmax))
    v = max(min(v, 255), 0)
    return palette[v]


def make_hash_sha256(o):
    "Thanks: https://stackoverflow.com/a/42151923"
    hasher = hashlib.sha256()
    hasher.update(repr(make_hashable(o)).encode())
    return hasher.hexdigest()


def make_hashable(o):
    "Thanks: https://stackoverflow.com/a/42151923"
    if isinstance(o, (tuple, list)):
        return tuple(sorted([make_hashable(e) for e in o]))

    if isinstance(o, dict):
        return tuple(sorted((k, make_hashable(v)) for k, v in o.items()))

    if isinstance(o, (set, frozenset)):
        return tuple(sorted(make_hashable(e) for e in o))
    return o


def UID(*args, length=7):
    chs = hashlib.sha512()
    for a in args:
        if isinstance(a, int) \
                or isinstance(a, float) \
                or isinstance(a, np.float32):
            chs.update(str(a).encode())
        elif isinstance(a, str):
            chs.update(str(a).lower().encode())
        elif isinstance(a, tuple):
            for x in a:
                chs.update(str(x).encode())
        elif isinstance(a, bytes):
            chs.update(a)
        elif isinstance(a, dict):
            for k, v in sorted(list(a.items())):
                chs.update(str(k).encode())
                chs.update(str(v).encode())
        elif isinstance(a, pd.Series):
            chs.update(a.values.tobytes())
        elif isinstance(a, pd.DataFrame):
            chs.update(a.values.tobytes())
        else:
            ft = str(type(a))
            raise Exception(f"Invalid type to generate UID {ft}")

    chs_hex = chs.hexdigest()
    return chs_hex[:length]


def exe(*cl, stdout=False):
    if stdout is True:
        P = sp.Popen(cl, stdout=sp.PIPE, text=True)
        out, err = P.communicate()
        assert P.returncode == 0
        return out
    else:
        P = sp.Popen(cl)
        P.communicate()
        return P.returncode


# thanks: https://stackoverflow.com/a/1094933
def sizeof_fmt(num, suffix='B'):
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


def diskcache(where='~/.cache/simple_disk_cache/',
              refresh=False, verbose=False):
    import inspect
    import pickle
    from pathlib import Path

    where = Path(where).expanduser()

    def dc(func):
        def wrapper(*args, **kwargs):
            serialized_args = [inspect.getsource(func)]
            serialized_args.extend(list(args))
            for k, v in sorted(kwargs.items()):
                serialized_args.extend([k, v])
            uid = UID(*serialized_args, length=64)
            cachedir = where / uid[:3] / uid[:6]
            cachefile = cachedir / uid

            if cachefile.exists() and not refresh:
                if verbose:
                    lg.info(f"return from cache - {uid}")

                with open(cachefile, 'rb') as F:
                    rv = pickle.load(F)
            else:
                if verbose:
                    lg.info(f"not in cache - {uid}")

                rv = func(*args, **kwargs)
                if not cachedir.exists():
                    cachedir.mkdir(parents=True)
                with open(cachefile, 'wb') as F:
                    pickle.dump(rv, F, protocol=4)

            return rv
        return wrapper
    return dc


def profiler(title, runs):
    import time

    def p2(func):
        def wrapper(*args, **kwargs):
            start = time.time()

            for x in range(runs):
                rv = func(*args, **kwargs)

            end = time.time()
            runtime = end - start
            rtitle = title[:30]
            print(f'$$$ {rtitle:40s}'
                  f' tt: {runtime:8.3g}'
                  f' no: {runs:4}'
                  f' 1: {(runtime/runs)*1000:8.4g} ms')
            return rv
        return wrapper
    return p2


@diskcache()
def query_pubmed(pubmed_id):
    from pymed import PubMed
    pubmed = PubMed(
        tool="Beehive", email="mark.fiers@vib.be")
    results = list(
        pubmed.query(f"{pubmed_id} [PMID]", max_results=100))
    r = results[0]

    rv = {}
    afs = "{firstname} {lastname}"
    if len(r.authors) == 1:
        rv['author'] = afs.format(**r.authors[0])
    else:
        rv['author'] = \
            (afs.format(**r.authors[0]) + ", " + afs.format(**r.authors[-1]))
    rv['title'] = r.title
    rv['doi'] = r.doi
    rv['year'] = r.publication_date.year
    return rv
