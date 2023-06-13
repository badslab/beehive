

import os
import logging
import duckdb
import pandas as pd
from typing import Optional


lg = logging.getLogger(__name__)

CONNECTION = None


def get_conn():
    global CONNECTION
    if CONNECTION is None:
        dbfile = os.environ['TERMITE_DB']
        CONNECTION = duckdb.connect(dbfile)
    return CONNECTION


def get_cursor():
    return get_conn().cursor()


def raw_sql(sql):
    conn = get_cursor()
    return conn.sql(sql).df()


def fuzzy_gene_search(experiment, datatype, gene):
    return raw_sql(f"""
        SELECT DISTINCT gene
          FROM expr
         WHERE experiment='{experiment}'
           AND datatype='{datatype}'
         ORDER BY jaro_winkler_similarity(gene, '{gene}') DESC
        LIMIT 25""")['gene']

    
def get_experiments():
    conn = get_cursor()
    sql = """SELECT * FROM help_experiments"""
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
    conn.sql("""DELET FROM experiments""")
    
    
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


def get_expr_obscat(experiment, datatype, gene, cat):
    sql = f"""
           SELECT expr.obs,
                  expr.value as expr,
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


def get_expr_obsnum(experiment, datatype, gene, num):
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
    print(sql)
    return raw_sql(sql)



def get_obs_category(experiment, cat):
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


def table_exists(conn, table):
    query = f"""
        SELECT EXISTS(
            SELECT 1 FROM information_schema.tables
             WHERE table_name = '{table}')"""
    result = conn.execute(query).fetchone()
    return result[0]


def table_count(table,
                conn: Optional[duckdb.DuckDBPyConnection] = None):
    if not table_exists(conn, table):
        return -1
    
    rv = conn.sql(f"SELECT count(*) FROM {table}").df().iloc[0,0]
    return rv
    

def helptables(conn: Optional[duckdb.DuckDBPyConnection] = None):

    if conn is None:
        conn = get_conn()

    lg.info("Create experiment help table")
    conn.sql("DROP TABLE IF EXISTS help_experiments")
    conn.sql("""
        CREATE TABLE help_experiments AS
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
        

def create_or_append(
        table: str,
        local_df: pd.DataFrame,
        conn: Optional[duckdb.DuckDBPyConnection] = None):
    
    #create a table - or if it exists - append
    if not table_exists(conn, table):
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
        t = table['table_name']
        c = table_count(table=t, conn=conn)
        rv[t] = c
    return rv
    

def drop_db(conn: Optional[duckdb.DuckDBPyConnection] = None):
    """Drop all tables."""
    if conn is None:
        conn = get_conn()
        
    for _, table in conn.sql('SHOW ALL TABLES').df().iterrows():
        print(f'drop {table}', table_count(conn, table))
        conn.sql(f'DROP TABLE IF EXISTS {table}')
    conn.execute("VACUUM")
    conn.execute("CHECKPOINT")

