"""Helper functions for beehive."""
import subprocess as sp
import hashlib

import beehive



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
                  default: str = None,
                  title: str = None,
                  curdoc=None,
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
    import random

    assert curdoc is not None

    param_name = name # to get & retrieve from the URL
    if title is None:
        title = name.capitalize()
        
    if widget == AutocompleteInput:
        # unique name - prevents browser based autocomplete
        # which messes up different datasets
        name = f"{name}.{random.randint(100000, 999999)}"

    js_onchange = CustomJS(
        code=f'insertUrlParam("{param_name}", cb_obj.value);')
    args = curdoc.session_context.request.arguments
    new_widget = widget(name=name, title=title, **kwargs)
    new_widget.value = getarg(args, param_name, default)
    new_widget.js_on_change("value", js_onchange)
    return new_widget


def getcolor(x, palette, vmax, vmin=None):
    if vmin is None:
        vmin = -vmax

    v = int(255 * (x + vmax) / (2 * vmax))
    v = max(min(v, 255), 0)
    return palette[v]


def UID(*args, length=7):
    chs = hashlib.sha512()
    for a in args:
        if isinstance(a, int):
            chs.update(str(a).encode())
        elif isinstance(a, str):
            chs.update(str(a).lower().encode())
        elif isinstance(a, bytes):
            chs.update(a)
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
            runtime = end-start
            rtitle = title[:30]
            print(f'$$$ {rtitle:40s}'
                  f' tt: {runtime:8.3g}'
                  f' no: {runs:4}'
                  f' 1: {(runtime/runs)*1000:8.4g} ms')
            return rv
        return wrapper
    return p2
