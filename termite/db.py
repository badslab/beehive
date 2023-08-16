

import logging
import os
from typing import Optional, List

import duckdb
import pandas as pd


lg = logging.getLogger(__name__)


CONNECTION: Optional[duckdb.DuckDBPyConnection] = None


def init_conn(dbfile: Optional[str] = None)  -> duckdb.DuckDBPyConnection:
    global CONNECTION
    
    if CONNECTION is not None:
        lg.debug("database already initialized?")
        exit(-1)
        
    if dbfile is None:
        dbfile = os.environ['TERMITE_DB']
    CONNECTION = duckdb.connect(dbfile)
    return CONNECTION


def get_conn() -> duckdb.DuckDBPyConnection:
    """ Return a global duckdb connection. One db per instance. """
    global CONNECTION
    if CONNECTION is None:
        lg.debug("database is not initialized")
        return init_conn()
    return CONNECTION



def get_cursor():
    """ Return a cursor to the connection."""
    return get_conn().cursor()


# @raw_sql_timer
def raw_sql(sql: str,
            conn: Optional[duckdb.DuckDBPyConnection] = None,
            ) -> pd.DataFrame:
    """Run SQL and return a Pandas Dataframe of the results."""
    
    if conn is None:
        conn = get_cursor()
        
    result = conn.sql(sql)
    
    if result is None:
        # empty dataframe in the case
        # the query returns nothing
        return pd.DataFrame([])
    else:
        return result.df()

def all_genes(exp_id: int) -> List[str]:

    md = raw_sql(
        f""" SELECT distinct gene
              FROM help_gene
             WHERE exp_id = {exp_id}
             ORDER BY gene""")

    return list(md['gene'])
    

# @st.cache_data() # (persist='disk')
def find_gene_candidates(
        exp_id: int) -> List[str]:

    md = raw_sql(
        f""" SELECT topgenes
              FROM experiment_md
             WHERE exp_id = {exp_id}
        """)

    if len(md) == 0:
        return []
    else:
        return md.iloc[0]['topgenes'].split()

    
# @st.cache_data() # (persist='disk')
def find_dimreds(
        exp_id: int) -> List[str]:

    md = raw_sql(
        f""" SELECT dimred
              FROM experiment_md
             WHERE exp_id = {exp_id}
        """)

    if len(md) == 0:
        return []
    else:
        return md.iloc[0]['dimred'].split()

    
# @st.cache_data
def fuzzy_gene_search(exp_id: int,
                      gene: str):
    return raw_sql(f"""
        SELECT DISTINCT gene
          FROM help_gene
         WHERE exp_id={exp_id}
         ORDER BY jaro_winkler_similarity(gene, '{gene}') DESC
        LIMIT 25""")['gene']

def get_experiments():
    return raw_sql('SELECT * FROM experiment_md')

def get_stats_experiment():
    conn = get_cursor()
    sql = """SELECT * FROM help_experiment"""
    return conn.sql(sql).df()


def get_obs_cat_names(exp_id: int):
    return raw_sql(f"""SELECT *
               FROM help_obs_cat
              WHERE exp_id={exp_id}
             LIMIT 100""")['name']


def get_obs_num_names(exp_id: int):
    return raw_sql(f"""SELECT *
               FROM help_obs_num
              WHERE exp_id={exp_id}
              LIMIT 100
    """)['name']


def forget(exp_id: int):
    conn = get_conn()
    for q in [
            f"DELETE FROM expr WHERE exp_id = {exp_id}",
            f"DELETE FROM experiment_md WHERE exp_id = {exp_id}",
            f"DELETE FROM obs_cat WHERE exp_id = {exp_id}",
            f"DELETE FROM obs_num WHERE exp_id = {exp_id}",
            f"DELETE FROM help_gene WHERE exp_id = {exp_id}",
            f"DELETE FROM help_experiment WHERE exp_id = {exp_id}",
            f"DELETE FROM help_obs_cat WHERE exp_id = {exp_id}",
            f"DELETE FROM help_obs_num WHERE exp_id = {exp_id}",
            ]:
        try:
            print(q)
            conn.sql(q)
        except duckdb.CatalogException:
            #assume the table does not exist?
            pass
    
    
# def get_expr_obscat_two(experiment, datatype, gene, cat1, cat2):
#     return raw_sql(f"""
#            SELECT expr.obs,
#                   expr.value as expr,
#                   occ1.value as cat1,
#                   occ2.value as cat2
#              FROM expr
#              JOIN obs_cat as occ1
#                ON (expr.experiment = occ1.experiment
#                    AND expr.obs = occ1.obs) 
#              JOIN obs_cat as occ2
#                ON (expr.experiment = occ2.experiment
#                    AND expr.obs = occ2.obs) 
#             WHERE expr.experiment='{experiment}'
#               AND expr.gene='{gene}'
#               AND expr.datatype='{datatype}'
#               AND occ1.name='{cat1}'
#               AND occ2.name='{cat2}'
#         """)


# def get_expr_obscat(
#         experiment: str,
#         datatype: str,
#         gene: str,
#         cat: str,
#         conn: Optional[duckdb.DuckDBPyConnection] = None):
    
#     sql = f"""
#            SELECT expr.obs,
#                   expr.value as num,
#                   occ1.value as cat,
#              FROM expr
#              JOIN obs_cat as occ1
#                ON (expr.experiment = occ1.experiment
#                    AND expr.obs = occ1.obs) 
#             WHERE expr.experiment='{experiment}'
#               AND expr.gene='{gene}'
#               AND expr.datatype='{datatype}'
#               AND occ1.name='{cat}'
#         """
#     return raw_sql(sql)


def get_obs_num(exp_id: int,
                name: str) -> pd.Series:

    return raw_sql(
        f""" SELECT obs, value
              FROM obs_num
             WHERE exp_id={exp_id}
               AND name='{name}'
        """).set_index('obs')['value']



def get_obs_cat(exp_id: int,
                name: str) -> pd.Series:

    return raw_sql(
        f""" SELECT obs, value
              FROM obs_cat
             WHERE exp_id={exp_id}
               AND name='{name}'
        """).set_index('obs')['value']


def get_expr_gene(exp_id: int,
                  gene: str) -> pd.Series:

    return raw_sql(
        f""" SELECT obs, value
              FROM expr
             WHERE exp_id={exp_id}
               AND gene='{gene}'
        """).set_index('obs')['value']


# def get_numvar_obscat(
#         exp_id: int,
#         numname: str,
#         catname: str,
#         conn: Optional[duckdb.DuckDBPyConnection] = None):
    
#     sql = f"""
#            SELECT num1.value as num,
#                   cat1.value as cat,
#              FROM obs_num as num1 
#              JOIN obs_cat as cat1
#                ON (num1.exp_id = cat1.exp_id
#                    AND num1.obs = cat1.obs) 
#             WHERE num1.exp_id={exp_id}
#               AND num1.name='{numname}'
#               AND cat1.name='{catname}'
#         """
#     return raw_sql(sql, conn=conn)


# def get_expr_obsnum(
#         exp_id: int,
#         gene: str,
#         num: str,
#         conn: Optional[duckdb.DuckDBPyConnection] = None):

#     if conn is None:
#         conn = get_conn()
        
#     sql = f"""
#            SELECT expr.obs,
#                   expr.value as expr,
#                   ocn1.value as num1,
#              FROM expr
#              JOIN obs_num as ocn1
#                ON (expr.experiment = ocn1.experiment
#                    AND expr.obs = ocn1.obs) 
#             WHERE expr.experiment='{experiment}'
#               AND expr.gene='{gene}'
#               AND expr.datatype='{datatype}'
#               AND ocn1.name='{num}'
#         """
#     return raw_sql(sql, conn=conn)



def get_obs_category(exp_id: int,
                     cat: str,
                     conn: Optional[duckdb.DuckDBPyConnection] = None):

    if conn is None:
        conn = get_cursor()
        
    sql = f"""SELECT obs_cat.value,
                     count(*) as count
               FROM obs_cat
              WHERE exp_id={exp_id}
                AND name='{cat}'
              GROUP BY obs_cat.value
              LIMIT 10
    """
    return conn.sql(sql).df()


def table_exists(table,
                 conn: Optional[duckdb.DuckDBPyConnection] = None) -> bool:
    
    if conn is None:
        conn = get_conn()
        
    rv = raw_sql(f"""
        SELECT EXISTS(
            SELECT 1 FROM information_schema.tables
             WHERE table_name = '{table}')""", conn=conn)
    return rv.iloc[0,0]



def table_count(table: str,
                conn: Optional[duckdb.DuckDBPyConnection] = None) -> int:

    if conn is None:
        conn = get_conn()
        
    if not table_exists(table, conn=conn):
        return -1
    
    rv = conn.sql(f"SELECT count(*) FROM {table}").df().iloc[0,0]
    return rv
    


def helptables(
        conn: Optional[duckdb.DuckDBPyConnection] = None,
        what: str = ""):

    if conn is None:
        conn = get_conn()

    if what and what in "help_dataset":
        lg.info("Create experiment help table")
        conn.sql("DROP TABLE IF EXISTS help_dataset")
        conn.sql("""
            CREATE TABLE help_dataset AS
              SELECT dataset_id,
                     count(*) as no_datapoints,
                     count(distinct obs) as no_cells,
                     count(distinct gene) as no_genes
                FROM expr
               GROUP BY dataset_id""")

    # if what and what in "help_obs_cat":
    #     lg.info("Create obs_cat help table")
    #     conn.sql("DROP TABLE IF EXISTS help_obs_cat")
    #     conn.sql("""
    #        CREATE TABLE help_obs_cat AS
    #            SELECT exp_id, name 
    #              FROM obs_cat
    #             GROUP BY exp_id, name """)

    # if what and what in "help_obs_num":
    #     lg.info("Create obs_num help table")
    #     conn.sql("DROP TABLE IF EXISTS help_obs_num")
    #     conn.sql("""
    #            CREATE TABLE help_obs_num AS
    #            SELECT DISTINCT exp_id, name
    #              FROM obs_num
    #             GROUP BY exp_id, name """)

    if what and what in "help_gene":
        lg.info("Create help_gene table")
        conn.sql("DROP TABLE IF EXISTS help_gene")
        conn.sql("""
           CREATE TABLE help_gene AS
               SELECT distinct dataset_id, gene,
                      SUM(value) AS sumval,
                      SUM(LEAST(value, 1)) / COUNT(value) AS fracnonzero
                 FROM expr
                GROUP BY dataset_id, gene
                ORDER BY dataset_id ASC, sumval DESC """)


def autoincrementor(table, field, value) -> int:
    try:
        result = raw_sql( f"""SELECT {field}_id
                                FROM {table}
                               WHERE {field} = '{value}' """)
        if len(result) > 0:
            return result.iloc[0,0]        
    except duckdb.CatalogException:
        # table does not exist?
        pass

    # no record
    try:
        max_id = raw_sql(
            f'''SELECT MAX({field}_id)
                FROM {table}''').iloc[0,0]
    except duckdb.CatalogException:
        # table does not exist?
        max_id = 0
    return max_id + 1
    

def create_or_append(
        table: str,
        local_df: pd.DataFrame | pd.Series | dict,
        conn: Optional[duckdb.DuckDBPyConnection] = None) -> None:

    if conn is None:
        conn = get_conn()

    if isinstance(local_df, dict):
        local_df = pd.DataFrame(pd.Series(local_df)).T
    elif isinstance(local_df, pd.Series):
        local_df = pd.DataFrame(local_df).T

    lg.debug(f"appending to {table} dataframe { local_df.shape }")
    local_df = local_df.sort_index(axis=1)
    
    #create a table - or if it exists - append
    if not table_exists(table):
        #lg("create & insert", table, local_df.shape)
        conn.sql(f"CREATE TABLE '{table}' AS SELECT * FROM local_df")
    else:
        #lg("append tbl", table, local_df.shape)
        sql = f"INSERT INTO '{table}' SELECT * FROM local_df"
        conn.sql(sql)

        
def all_table_count(conn: Optional[duckdb.DuckDBPyConnection] = None):
    if conn is None:
        conn = get_conn()
    
    rv = {}
    for _, table in conn.sql('SHOW ALL TABLES').df().iterrows():
        t = table.loc['name']
        c = table_count(table=t)
        rv[t] = c
    return rv

    

def drop_db(conn: Optional[duckdb.DuckDBPyConnection] = None):
    """Drop all tables."""
    if conn is None:
        conn = get_conn()
        
    for _, table in conn.sql('SHOW ALL TABLES').df().iterrows():
        norec = table_count(table.loc["name"])
        print(f'Drop {table} (no rec: {norec}')
        conn.sql(f'DROP TABLE IF EXISTS {table}')
    conn.execute("VACUUM")
    conn.execute("CHECKPOINT")

