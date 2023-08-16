
from os import name
from termite.vis import util
from termite.util import cache


def get_dataset_stats(ds_id):
    return util.SqlOneRecord(
        sql=f"""
            SELECT * FROM help_dataset where dataset_id = {ds_id}
        """ )

def get_obs_data_exp(exp_id):
    return util.execute_sql(
        sql=f"""
        SELECT *
          FROM help_obs
         WHERE exp_id={exp_id}
        """)

def get_categ_data_exp(exp_id):
    return util.execute_sql(
        sql=f"""
        SELECT *
          FROM help_obs
         WHERE exp_id={exp_id}
           AND dtype='categorical'
        """)


@cache.memoize(typed=True, expire=3600)
def get_obs_cat(exp_id, name):
    rv = util.execute_sql(
        sql=f"""
        SELECT obs, value
          FROM obs_cat
         WHERE exp_id={exp_id}
           AND name='{name}'
        """).set_index('obs')['value']
    rv.name = name
    return rv


@cache.memoize(typed=True, expire=3600)
def get_obs_num(exp_id, name):
    rv = util.execute_sql(
        sql=f"""
        SELECT obs, value
          FROM obs_num
         WHERE exp_id={exp_id}
           AND name='{name}'
        """).set_index('obs')['value']
    rv.name = name
    return rv


@cache.memoize(typed=True, expire=3600)
def get_expr(ds_id, gene):
    rv = util.execute_sql(
        f"""SELECT obs, value
              FROM expr
             WHERE dataset_id = {ds_id}
               AND gene = '{gene}' """).set_index('obs')['value']
    rv.name = gene
    return rv


@cache.memoize(typed=True, expire=3600)
def get_topgenes(ds_id):
    return util.execute_sql(
        sql=f"""
            SELECT * FROM help_gene where dataset_id = {ds_id}
             ORDER BY sumval DESC
             LIMIT 20
        """ )


# def get_dataset_stats(ds_id):
#     db = get_conn()
#     return util.SqlOneRecord(
#         db=db,
#         sql=f"""
#             SELECT * FROM help_dataset where dataset_id = {ds_id}
#         """ )


# def get_obs_data_exp(exp_id):
#     db = get_conn()
#     return util.execute_sql(
#         db=db,
#         sql=f"""
#         SELECT *
#           FROM help_obs
#          WHERE exp_id={exp_id}
#         """)

# def get_categ_data_exp(exp_id):
#     db = get_conn()
#     return util.execute_sql(
#         db=db,
#         sql=f"""
#         SELECT *
#           FROM help_obs
#          WHERE exp_id={exp_id}
#            AND dtype='categorical'
#         """)


# def get_categ(exp_id, name):
#     return util.execute_sql(
#         sql=f"""
#         SELECT *
#           FROM obs_cat
#          WHERE exp_id={exp_id}
#            AND name='{name}'
#         """)
    
    
# def get_topgenes(ds_id):
#     db = get_conn()
#     return util.execute_sql(
#         db=db,
#         sql=f"""
#             SELECT * FROM help_gene where dataset_id = {ds_id}
#              ORDER BY sumval DESC
#              LIMIT 20
#         """ )


# def get_dataset_stats(ds_id):
#     db = get_conn()
#     return util.SqlOneRecord(
#         db=db,
#         sql=f"""
#             SELECT * FROM help_dataset where dataset_id = {ds_id}
#         """ )


# def get_obs_data_exp(exp_id):
#     db = get_conn()
#     return util.execute_sql(
#         db=db,
#         sql=f"""
#         SELECT *
#           FROM help_obs
#          WHERE exp_id={exp_id}
#         """)

# def get_categ_data_exp(exp_id):
#     db = get_conn()
#     return util.execute_sql(
#         db=db,
#         sql=f"""
#         SELECT *
#           FROM help_obs
#          WHERE exp_id={exp_id}
#            AND dtype='categorical'
#         """)


# def get_categ(exp_id, name):
#     return util.execute_sql(
#         sql=f"""
#         SELECT *
#           FROM obs_cat
#          WHERE exp_id={exp_id}
#            AND name='{name}'
#         """)
    
    
# def get_topgenes(ds_id):
#     db = get_conn()
#     return util.execute_sql(
#         db=db,
#         sql=f"""
#             SELECT * FROM help_gene where dataset_id = {ds_id}
#              ORDER BY sumval DESC
#              LIMIT 20
#         """ )
