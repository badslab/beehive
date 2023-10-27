

import gzip
import logging
from functools import lru_cache, partial
from pathlib import Path

lg = logging.getLogger(__name__)


def simple_disk_cache(cache_path: Path,
                      cache_name: str,
                      force: bool = False,
                      lru:int = 16):
    """Very simple decorator to cache functions output to disk

    """

    import copy
    import inspect
    import pickle
    import re
    cache_name = cache_name.strip('/')
    if not cache_path.exists():
        cache_path.mkdir(parents=True)

    def decorator(func):
        """The decorator funtion."""
        def _md5(x):
            from hashlib import md5
            m = md5()
            m.update(str(x).encode())
            return m.hexdigest()

        #lru_cache to prevent rehashing the same object too much
        # not sure if this really helps - as I expect lru_cache
        # hashes as well...
        @lru_cache(256)
        def _umd5(x):
            from hashlib import md5
            m = md5()
            m.update(str(x).upper().encode())
            return m.hexdigest()

        # embedded functions
        varfunc = {
            'md5': _md5,
            'umd5': _umd5,
            'hash': hash,
            'upper': str.upper,
            'lower': str.lower }

        _cache_path = cache_path / func.__name__
        if not _cache_path.exists():
            _cache_path.mkdir(parents=True)

        # if you want the cache to depend on the
        # source code as well
        src_hash = abs(hash(inspect.getsource(func)))
        #@lru_cache(maxsize=lru)
        def wrapper(*args, **kwargs):
            """And the wrapped function."""

            finc = inspect.getfullargspec(func)

            if finc.defaults is None:
                defkwargs = {}
            else:
                defkwargs = dict(zip(
                    finc.args[-len(finc.defaults):], finc.defaults))

            full_kwargs = copy.copy(kwargs)
            for k, v in zip(inspect.getfullargspec(func).args, args):
                full_kwargs[k] = v
            for k, v in defkwargs.items():
                if not k in full_kwargs:
                    full_kwargs[k] = v

            #apply possible functions to variable name
            rex = re.compile(r"\{([^}]+?)\|([A-Za-z0-9\|]+)\}")
            for var, fname in rex.findall(cache_name):
                if not var in full_kwargs:
                    raise KeyError(f"did not find var: {var}")
                newval = full_kwargs[var]
                for f in fname.split('|'):
                    assert f in varfunc
                    newval = varfunc[f](newval)
                #print(var, fname, newval)
                full_kwargs[f"{var}|{fname}"] = newval

            #for m in fm.iter(fm)
            _cache_file = _cache_path / \
                cache_name.format(CODE=src_hash, **full_kwargs)

            print(_cache_file)

            if not _cache_file.parent == _cache_path:
                if not _cache_file.parent.exists():
                    _cache_file.parent.mkdir(parents=True)
            if not force and _cache_file.exists():
                with gzip.open(_cache_file, mode='rb') as F:
                    result = pickle.load(F)
            else:
                result = func(*args, **kwargs)
                with gzip.open(_cache_file, compresslevel=3, mode='wb') as F:
                    pickle.dump(result, F)
            return result
        return wrapper
    return decorator

sdc_cache_path = Path('~/.cache/termite').expanduser().resolve()
sdc = partial(simple_disk_cache, sdc_cache_path)


@sdc("{pubmed_id}_{CODE}")
def query_pubmed(pubmed_id):
    """Get pubmed data."""
    from pymed import PubMed
    pubmed = PubMed(
        tool="Termite", email="mark.fiers@vib.be")
    results = list(
        pubmed.query(f"{pubmed_id} [PMID]", max_results=100))
    r = results[0]

    rv = {}
    afs = "{firstname} {lastname}"
    if len(r.authors) == 1:
        rv['author'] = afs.format(**r.authors[0])
    elif len(r.authors) == 2:
        rv['author'] = "{} and {}".format(
            afs.format(**r.authors[0]), afs.format(**r.authors[0]))
    else:
        rv['author'] = \
            (afs.format(**r.authors[0])
             + " ... "
             + afs.format(**r.authors[-1]))
    rv['title'] = r.title
    rv['abstract'] = r.abstract
    rv['doi'] = r.doi
    rv['year'] = r.publication_date.year
    return rv




# dbfile = Path(os.environ['TERMITE_DB'])
# cache_folder = dbfile.parent / f"{dbfile.name}_diskcache"
# if not cache_folder.exists():
#     cache_folder.mkdir()
# lg.warning(f"caching to {cache_folder}")
# cache = FanoutCache(directory=str(cache_folder), expire=60*60)


# def prep_sql(sql, data):
#     sql = sql.format(**data)
#     return sql

# def execute_sql(item, dip):
#     sql = item['sql']
#     sql = prep_sql(sql, dip.data)
#     lg.debug(f"executing sql:\n{sql}\n--------------")
#     df = dip.conn.execute(sql).df()
#     lg.debug(f"result: {len(df)} records")
#     return df


# def raw_sql_timer(func):
#     from sql_formatter.core import format_sql
#     @wraps(func)
#     def timeit_wrapper(*args, **kwargs):
#         start_time = time.perf_counter()
#         print('#### START SQL ################################################')
#         print(format_sql(args[0]))
#         result = func(*args, **kwargs)
#         end_time = time.perf_counter()
#         total_time = end_time - start_time
#         print(f'######### READY {total_time:12.4f} sec / {result.shape}')
#         return result
#     return timeit_wrapper
