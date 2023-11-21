

from typing import Optional, TYPE_CHECKING

from .db import CHDB


if TYPE_CHECKING:
    import pandas as pd


class API():
    def __init__(self,
                 dbfile: Optional[str] = None,
                 read_only: bool = True) -> None:
        self.dbfile = dbfile
        self.db = CHDB(dbfile, read_only)

    def datasets(self, query: str) -> "pd.DataFrame":
        sql = """
            SELECT DISTINCT *
              FROM experiment_md """
        if query:
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
