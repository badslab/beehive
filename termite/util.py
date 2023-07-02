
import hashlib
from functools import wraps
import logging
import time


lg = logging.getLogger(__name__)


def raw_sql_timer(func):
    from sql_formatter.core import format_sql
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        print('#### START SQL ################################################')
        print(format_sql(args[0]))
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'######### READY {total_time:12.4f} sec / {result.shape}')
        return result
    return timeit_wrapper


def UID(*args, length=7):
    chs = hashlib.sha512()
    import numpy as np
    import pandas as pd

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
