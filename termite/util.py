
import hashlib
from functools import wraps
import logging
from pathlib import Path
import os
import time


import duckdb


from diskcache import FanoutCache

lg = logging.getLogger(__name__)

dbfile = Path(os.environ['TERMITE_DB'])
cache_folder = dbfile.parent / f"{dbfile.name}_diskcache"
if not cache_folder.exists():
    cache_folder.mkdir()
lg.warning(f"caching to {cache_folder}")
cache = FanoutCache(directory=str(cache_folder), expire=60*60)


def prep_sql(sql, data):
    sql = sql.format(**data)
    return sql

def execute_sql(item, dip):
    sql = item['sql']
    sql = prep_sql(sql, dip.data)
    lg.debug(f"executing sql:\n{sql}\n--------------")
    df = dip.conn.execute(sql).df()
    lg.debug(f"result: {len(df)} records")
    return df


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


def query_pubmed(pubmed_id):
    from pymed import PubMed
    pubmed = PubMed(
        tool="Beehive", email="mark.fiers@vib.be")
    results = list(
        pubmed.query(f"{pubmed_id} [PMID]", max_results=100))
    r = results[0]

    print(r)
    rv = {}
    afs = "{firstname} {lastname}"
    if len(r.authors) == 1:
        rv['author'] = afs.format(**r.authors[0])
    else:
        rv['author'] = \
            (afs.format(**r.authors[0])
             + ", "
             + afs.format(**r.authors[-1]))
    rv['title'] = r.title
    rv['abstract'] = r.abstract
    rv['doi'] = r.doi
    rv['year'] = r.publication_date.year
    return rv
