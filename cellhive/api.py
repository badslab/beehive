

from typing import TYPE_CHECKING, List, Optional

from .db import CHDB

if TYPE_CHECKING:
    import pandas as pd


class API():
    def __init__(self,
                 dbfile: Optional[str] = None,
                 read_only: bool = True) -> None:
        self.dbfile = dbfile
        self.db = CHDB(dbfile=dbfile,
                       read_only=read_only)

    def sql(self, sql: str):
        """
        This method allows you to pass a SQL query (in the form of a string)
        to be executed against the underlying database. The results are
        returned by the `db` object.

        Parameters:
        - sql (str): The SQL query to execute.

        Returns:
        - A pandas dataframe with the results from executing the SQL
          statement on the database.
        """
        return self.db.sql(sql)


    def obsm_names(self, exp_id: int):
        sql = f"""
            SELECT DISTINCT name
            FROM obs_num
            WHERE exp_id = {exp_id}
              AND name LIKE '%/0'
        """
        rv = self.db.sql(sql)
        rv = [x.replace('/0', '') for x in rv['name']]
        return rv


    def obsm(self, exp_id: int, obsm_name: str):
        sql = f"""
            SELECT *
            FROM obs_num
            WHERE exp_id = {exp_id}
              AND name LIKE '{obsm_name}/%'
        """
        rv = self.db.sql(sql)
        rv = rv.pivot(index='cell', columns='name', values='value')
        rv = rv.sort_index(axis=1, ascending=True)
        return rv


    def obs_num(self, exp_id: int, obs_names: List[str] | str):

        if isinstance(obs_names, str):
            obs_names = [obs_names]

        inlist = ( "('"
                   + "', '".join(obs_names)
                   + "')" )
        sql = f"""
            SELECT *
            FROM obs_num
            WHERE exp_id = {exp_id}
              AND name IN {inlist} """
        rv = self.db.sql(sql)
        rv = rv.pivot(index='cell', columns='name', values='value')
        rv = rv.sort_index(axis=1, ascending=True)
        return rv

    def obs_names_num(self, exp_id: int):
        sql = f"""
            SELECT DISTINCT name
            FROM obs_num
            WHERE exp_id = {exp_id} """
        rv = self.db.sql(sql)
        return rv



    def obs_cat(self, exp_id: int, obs_names: List[str] | str):

        if isinstance(obs_names, str):
            obs_names = [obs_names]

        inlist = ( "('"
                   + "', '".join(obs_names)
                   + "')" )
        sql = f"""
            SELECT *
            FROM obs_cat
            WHERE exp_id = {exp_id}
              AND name IN {inlist} """
        rv = self.db.sql(sql)
        rv = rv.pivot(index='cell', columns='name', values='value')
        rv = rv.sort_index(axis=1, ascending=True)
        return rv


    def obs_names_cat(self, exp_id: int):
        sql = f"""
            SELECT DISTINCT name
            FROM obs_cat
            WHERE exp_id = {exp_id} """
        rv = self.db.sql(sql)
        return rv


    def datasets(self, query: str | None = None) -> "pd.DataFrame":
        sql = """
            SELECT DISTINCT *
              FROM experiment_md """
        if query is not None:
            sql += f"""
             WHERE ( title ILIKE '%{query}%'
                     OR author ILIKE  '%{query}%'
                     OR abstract ILIKE '%{query}%' )
            """
        rv = self.db.sql(sql)
        return rv


    def gene(self, dataset_id: int, gene: str) -> "pd.Series":
        sql = f"""
            SELECT obs, value
              FROM expr
            WHERE dataset_id = {dataset_id}
        """
        rv = self.db.sql(sql).set_index('obs')['value']
        rv.name = gene
        return rv
