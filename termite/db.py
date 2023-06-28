

import os
import logging
import duckdb
import pandas as pd
from typing import Optional, List
import streamlit as st


lg = logging.getLogger(__name__)

CONNECTION: Optional[duckdb.DuckDBPyConnection] = None


def get_conn() -> duckdb.DuckDBPyConnection:
    """ Return a global duckdb connection. One db per instance. """
    global CONNECTION
    if CONNECTION is None:
        dbfile = os.environ['TERMITE_DB']
        CONNECTION = duckdb.connect(dbfile)
    return CONNECTION


def get_cursor():
    """ Return a cursor to the connection."""
    return get_conn().cursor()


def raw_sql(sql: str,
            conn: Optional[duckdb.DuckDBPyConnection] = None) -> pd.DataFrame:
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

 
@st.cache_data() # (persist='disk')
def find_gene_candidates(
        experiment: str,
        datatype: str) -> List[str]:


    favorites = '''APOE APP PSEN1 PSEN2 CLU GFAP PLP1
                   BIN1 SORL1 ABCA7 CR1 CD33 MS4A
                   PICALM EPHA1 INPPD5 MEF2C ADAM10
                   AKAP4 MEG3'''.lower().split()

    never = '''MALAT1 MT-CO1'''.lower().split()
    
    favlist = "('" + "','".join(favorites) + "')"
    neverlist = "('" + "','".join(never) + "')"

    c1 = raw_sql(f'''
        SELECT gene, QUANTILE_CONT(value, 0.99) as medexpr FROM expr
         WHERE experiment = '{experiment}'
          AND datatype = '{datatype}'
          AND lower(gene) IN {favlist}
        GROUP BY gene
       HAVING (medexpr > 0)
        ORDER BY medexpr DESC
        LIMIT 10
    ''')
    
    c2 = raw_sql(f'''
        SELECT gene, QUANTILE_CONT(value, 0.99) as medexpr
          FROM expr
         WHERE experiment = '{experiment}'
          AND datatype = '{datatype}'
          AND lower(gene) NOT IN {neverlist}
          AND NOT starts_with(lower(gene), 'mt-')
     
        GROUP BY gene
        HAVING (medexpr > 0)
        ORDER BY medexpr DESC
        LIMIT 10
    ''')

    
    c1['medexpr'] = c1['medexpr'] * 10
    
    rv = list(pd.concat([c1, c2]).sort_values(by='medexpr')['gene'])
    lg.info(f'Candidate genes: {rv[:5]}')
    
    return rv

    
@st.cache_data
def fuzzy_gene_search(experiment, datatype, gene):
    return raw_sql(f"""
        SELECT DISTINCT gene
          FROM help_gene
         WHERE experiment='{experiment}'
           AND datatype='{datatype}'
         ORDER BY jaro_winkler_similarity(gene, '{gene}') DESC
        LIMIT 25""")['gene']

    
def get_exp_datatypes():
    conn = get_cursor()
    sql = """SELECT * FROM help_datatype"""
    return conn.sql(sql).df()


def get_obs_cat_names(experiment):
    return raw_sql(f"""SELECT *
               FROM help_obs_cat
              WHERE experiment='{experiment}'
             LIMIT 100""")['name']


def get_obs_num_names(experiment):
    return raw_sql(f"""SELECT *
               FROM help_obs_num
              WHERE experiment='{experiment}'
             LIMIT 100""")['name']


def forget(experiment):
    conn = get_conn()
    for q in [
            f"DELETE FROM expr WHERE experiment = '{experiment}'",
            f"DELETE FROM experiment_md WHERE experiment = '{experiment}'",
            f"DELETE FROM obs_cat WHERE experiment = '{experiment}'",
            f"DELETE FROM obs_num WHERE experiment = '{experiment}'",
            f"DELETE FROM help_datatypes WHERE experiment = '{experiment}'",
            f"DELETE FROM help_obs_cat WHERE experiment = '{experiment}'",
            f"DELETE FROM help_obs_num WHERE experiment = '{experiment}'",
            ]:
        try:
            conn.sql(q)
        except duckdb.CatalogException:
            #assume the table does not exist?
            pass
    
    
def get_expr_obscat_two(experiment, datatype, gene, cat1, cat2):
    return raw_sql(f"""
           SELECT expr.obs,
                  expr.value as expr,
                  occ1.value as cat1,
                  occ2.value as cat2
             FROM expr
             JOIN obs_cat as occ1
               ON (expr.experiment = occ1.experiment
                   AND expr.obs = occ1.obs) 
             JOIN obs_cat as occ2
               ON (expr.experiment = occ2.experiment
                   AND expr.obs = occ2.obs) 
            WHERE expr.experiment='{experiment}'
              AND expr.gene='{gene}'
              AND expr.datatype='{datatype}'
              AND occ1.name='{cat1}'
              AND occ2.name='{cat2}'
        """)


def get_expr_obscat(
        experiment: str,
        datatype: str,
        gene: str,
        cat: str,
        conn: Optional[duckdb.DuckDBPyConnection] = None):
    
    sql = f"""
           SELECT expr.obs,
                  expr.value as num,
                  occ1.value as cat,
             FROM expr
             JOIN obs_cat as occ1
               ON (expr.experiment = occ1.experiment
                   AND expr.obs = occ1.obs) 
            WHERE expr.experiment='{experiment}'
              AND expr.gene='{gene}'
              AND expr.datatype='{datatype}'
              AND occ1.name='{cat}'
        """
    return raw_sql(sql)


def get_obs_num(experiment: str,
                name: str) -> pd.Series:

    return raw_sql(
        f""" SELECT obs, value
              FROM obs_num
             WHERE experiment='{experiment}'
               AND name='{name}'
        """).set_index('obs')['value']



def get_obs_cat(experiment: str,
                name: str) -> pd.Series:

    return raw_sql(
        f""" SELECT obs, value
              FROM obs_cat
             WHERE experiment='{experiment}'
               AND name='{name}'
        """).set_index('obs')['value']


def get_expr_gene(experiment: str,
                  datatype: str,
                  gene: str) -> pd.Series:

    return raw_sql(
        f""" SELECT obs, value
              FROM expr
             WHERE experiment='{experiment}'
               AND datatype='{datatype}'
               AND gene='{gene}'
        """).set_index('obs')['value']


def get_numvar_obscat(
        experiment: str,
        numname: str,
        catname: str,
        conn: Optional[duckdb.DuckDBPyConnection] = None):
    
    sql = f"""
           SELECT num1.value as num,
                  cat1.value as cat,
             FROM obs_num as num1 
             JOIN obs_cat as cat1
               ON (num1.experiment = cat1.experiment
                   AND num1.obs = cat1.obs) 
            WHERE num1.experiment='{experiment}'
              AND num1.name='{numname}'
              AND cat1.name='{catname}'
        """
    return raw_sql(sql, conn=conn)


def get_expr_obsnum(
        experiment: str,
        datatype: str,
        gene: str,
        num: str,
        conn: Optional[duckdb.DuckDBPyConnection] = None):

    if conn is None:
        conn = get_conn()
        
    sql = f"""
           SELECT expr.obs,
                  expr.value as expr,
                  ocn1.value as num1,
             FROM expr
             JOIN obs_num as ocn1
               ON (expr.experiment = ocn1.experiment
                   AND expr.obs = ocn1.obs) 
            WHERE expr.experiment='{experiment}'
              AND expr.gene='{gene}'
              AND expr.datatype='{datatype}'
              AND ocn1.name='{num}'
        """
    return raw_sql(sql, conn=conn)



def get_obs_category(experiment, cat,
                     conn: Optional[duckdb.DuckDBPyConnection] = None):

    if conn is None:
        conn = get_cursor()
        
    sql = f"""SELECT obs_cat.value,
                     count(*) as count
               FROM obs_cat
              WHERE experiment='{experiment}'
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
    


def helptables(conn: Optional[duckdb.DuckDBPyConnection] = None):

    if conn is None:
        conn = get_conn()

    lg.info("Create datatype help table")
    conn.sql("DROP TABLE IF EXISTS help_datatype")
    conn.sql("""
        CREATE TABLE help_datatype AS
          SELECT experiment, datatype, 
                 count(*) as no_datapoints,
                 count(distinct obs) as no_cells,
                 count(distinct gene) as no_genes
            FROM expr
           GROUP BY experiment, datatype""")

    lg.info("Create obs_cat help table")
    conn.sql("DROP TABLE IF EXISTS help_obs_cat")
    conn.sql("""
       CREATE TABLE help_obs_cat AS
           SELECT experiment, name 
             FROM obs_cat
            GROUP BY experiment, name """)

    lg.info("Create obs_num help table")
    conn.sql("DROP TABLE IF EXISTS help_obs_num")
    conn.sql("""
       CREATE TABLE help_obs_num AS
           SELECT experiment, name 
             FROM obs_num
            GROUP BY experiment, name """)
    
    lg.info("Create obs_num help table")
    conn.sql("DROP TABLE IF EXISTS help_gene")
    conn.sql("""
       CREATE TABLE help_gene AS
           SELECT distinct experiment, datatype, gene
             FROM expr """)


def create_or_append(
        table: str,
        local_df: pd.DataFrame,
        conn: Optional[duckdb.DuckDBPyConnection] = None) -> None:

    if conn is None:
        conn = get_conn()

        lg.debug(f"appending to {table} dataframe { local_df.shape }")
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
        t = table['name']
        c = table_count(table=t)
        rv[t] = c
    return rv
    

def drop_db(conn: Optional[duckdb.DuckDBPyConnection] = None):
    """Drop all tables."""
    if conn is None:
        conn = get_conn()
        
    for _, table in conn.sql('SHOW ALL TABLES').df().iterrows():
        print(f'Drop {table} (no rec: {table_count(table)}')
        conn.sql(f'DROP TABLE IF EXISTS {table}')
    conn.execute("VACUUM")
    conn.execute("CHECKPOINT")

