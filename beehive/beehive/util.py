"""Helper functions for beehive."""
import subprocess as sp
import hashlib
import time
from typing import List
import pandas as pd
import beehive
import numpy as np

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
        print(f'Function run {1000*(t2-t1):.2f}ms.')
        return res
    return wrapper


def get_datadir(name):
    """Return a bokeh view's data folder."""
    return beehive.DATADIR / name


def getarg(args, name, default=None, dtype=str):
    """Bokeh helper to get arguments from the URL."""
    if name in args:
        val = args[name][0].decode()
        return dtype(val)
    else:
        return default


def create_widget(name: str,
                  widget,
                  default=None,
                  title: str = None,
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
    from bokeh.models.callbacks import CustomJS
    from bokeh.models.widgets.inputs import AutocompleteInput
    from bokeh.models import RadioGroup
    import random

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
        new_widget = widget(name=name,title = title, **kwargs)


    if (value_type == list) \
            and not(type(getarg(args, param_name, default)) == list):
        new_widget.value = getarg(args, param_name, default).split(",")
    elif (value_type == int) or (value_type == float):
        new_widget.value = getarg(args, param_name, default, dtype=value_type)
    else:
        new_widget.value = getarg(args, param_name, default)
    if update_url:
        new_widget.js_on_change("value", js_onchange_value)
    return new_widget


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
        if isinstance(a, int) or isinstance(a,float) or isinstance(a,np.float32):
            chs.update(str(a).encode())
        elif isinstance(a, str):
            chs.update(str(a).lower().encode())
        elif isinstance(a, bytes):
            chs.update(a)
        elif isinstance(a,pd.Series):
            chs.update(str(a).encode())
        else:
            raise Exception("Invalid type to generate UID")
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


def diskcache(where='~/.cache/simple_disk_cache/', refresh=False):
    from pathlib import Path
    import inspect
    import pickle

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
                with open(cachefile, 'rb') as F:
                    rv = pickle.load(F)
            else:
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
