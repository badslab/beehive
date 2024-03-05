"""Database related code."""

import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Union

from typing_extensions import LiteralString

if TYPE_CHECKING:
    import pandas as pd
    from anndata import AnnData


lg = logging.getLogger()


class CHDB:
    """
    One class wraps many db functions.

    Prevents from having to pass the connection around all the time.
    """

    EXPERIMENT_MD_KEYS = """abstract author dataset_id
        dataset doi experiment experiment_id full_experiment genes
        layer_name organism pubmed study study_id title version year""".split()

    def __init__(self,
                 dbfile: Optional[str] = None,
                 read_only: bool = False) -> None:
        """Construct the chdb object."""

        import duckdb

        if dbfile is None:
            if 'CELLHIVE_DB' in os.environ:
                dbfile = os.environ['CELLHIVE_DB']
                if not Path(dbfile).exists():
                    raise FileNotFoundError("Duckdb not found")
            else:
                dbfile = 'cellhive.duckdb'

        self.dbfile = dbfile
        self.read_only = read_only
        self.conn = duckdb.connect(self.dbfile, read_only=read_only)

    def status(self) -> Dict[str, Any]:
        """Return a few db statistics."""
        rv = {}
        rv['dbsize'] = os.path.getsize(self.dbfile)
        for table, tcount in self.all_table_count().items():
            rv[table] = tcount
        return rv


    def table_exists(self, table: str) -> bool:
        """Check if a table exists."""
        rv = self.sql(f"""
            SELECT EXISTS(
                SELECT 1 FROM information_schema.tables
                 WHERE table_name = '{table}')""")

        return bool(rv.iloc[0,0])


    def all_table_count(self) -> Dict[str, int]:
        rv = {}
        for _, table_info in self.sql('SHOW ALL TABLES').iterrows():
            table = table_info['name']
            rv[table] = self.table_count(table)
        return rv

    def table_count(self, table: str) -> int:
        #if not table_exists(table, conn=conn):
        #    return -1
        rv = self.sql(f"SELECT count(*) as cnt FROM {table}")
        return rv['cnt'].iloc[0]


    def sql(self, sql: str) -> "pd.DataFrame":
        """Run SQL and return a Pandas Dataframe of the results."""

        import pandas as pd
        result = self.conn.sql(sql)

        if result is None:
            # empty dataframe in the case
            # the query returns nothing
            return pd.DataFrame([])
        return result.df()

    def get_id(self,
               table: str,
               field: str,
               value: Any
               ) -> int:
        """
        This function auto-increments an ID in a specific field of a
        database table.

        Parameters:
        table (str): The name of the table in the database to interact
            with.
        field (str): The specific field within the table to increment.
        value (str): The value to search for within the chosen field.

        The function first tries to select the ID related to the
        provided value in the specified field. If it exists, the
        function will return this ID.

        If the table doesn't exist or there's no record corresponding to
        the provided value, the function then tries to retrieve the
        maximum ID in the mentioned field and returns this ID
        incremented by 1. If the table still does not exist, it will
        return 1 (0 + 1) as the new ID to use.

        Returns:
        int: The ID to use next in the field.
        """

        import duckdb
        try:
            sql = f"""SELECT DISTINCT {field}_id
                        FROM {table}
                       WHERE {field} = '{value}' """
            result = self.sql(sql)

            if len(result) > 0:
                return result.iloc[0,0]

        except duckdb.CatalogException:
            # table does not exist?
            pass

        # no record
        try:
            max_id = self.sql(
                f'''SELECT MAX({field}_id)
                    FROM {table}''').iloc[0,0]
        except duckdb.CatalogException:
            # table does not exist?
            max_id = 0
        return max_id + 1


    def import_count_table(self,
                           dataset_id: int,
                           adata: "AnnData",
                           layer: str,
                           chunksize: int = 50000,
                           ) -> None:
        """Import an adata count matrix."""
        conn = self.conn

        lg.info("Start storing expression matrix")
        lg.info(f"Processing layer {layer}")

        if layer == 'X':
            x = adata.to_df()
        elif layer == 'RAW':
            x = adata.raw.to_adata().to_df()
        else:
            x = adata.to_df(layer=layer)

        # chunk this - using too much memory:
        x.index.name = 'obs'
        x.columns.name = 'gene'

        #if normalize == 'logrpm':
        #    x = np.log1p(10000 * x.copy().divide(x.sum(1), axis=0))


        #remove old data
        lg.info("remove old data")
        if self.table_exists('expr'):
            sql = f"""
                DELETE FROM expr
                 WHERE dataset_id={dataset_id}"""
            conn.sql(sql)

        lg.info("start expression data upload")
        for ic in range(0, x.shape[0], chunksize):
            chunk = x.iloc[ic:ic+chunksize,:]
            # chunk_rnk = rnk.iloc[ic:ic+chunksize,:]

            melted = chunk.reset_index().melt(id_vars='obs')
            melted['value'] = melted['value'].astype(float)  # ensure!
            melted['dataset_id'] = dataset_id

            # melted_rnk = chunk.reset_index().melt(id_vars='obs', value_name='rank')
            # melted_rnk['rank'] = melted_rnk['rank'].astype(float)  # ensure!

            # to be sure!
            # assert melted[['obs', 'gene']].equals(melted_rnk[['obs', 'gene']])
            # melted['rank'] = melted_rnk['rank']

            lg.info(f"chunk {ic}/{x.shape[0]} - melt {melted.shape[0]:_d}")
            self.create_or_append('expr', melted)

        # ensure index
        # print(db.raw_sql('create index idx_expr_eg on expr (exp_id, gene)'))


    def _check_incoming_df(self,
                           local_df: Union["pd.DataFrame", "pd.Series", dict],
                           ) -> "pd.DataFrame":

        """Check incoming df, force to dataframe."""

        import pandas as pd
        if isinstance(local_df, dict):
            local_df = pd.DataFrame(pd.Series(local_df)).T
        elif isinstance(local_df, pd.Series):
            local_df = pd.DataFrame(local_df).T

        local_df = local_df.sort_index(axis=1)
        return local_df


    def store_obscol(self,
                     col: "pd.Series",
                     name: str,
                     dtype: str,
                     exp_id: str):
        """Store an obs/obsm columns."""
        import pandas as pd

        col.index.name = 'cell'
        col.name = 'value'

        d = pd.DataFrame(col).reset_index()
        d['name'] = name
        d['exp_id'] = exp_id

        if dtype == 'cat':
            d['value'] = d['value'].astype(str)
            table = 'obs_cat'

        elif dtype in ['int', 'float', 'dimred']:
            d['value'] = d['value'].astype(float)
            table = 'obs_num'
        else:
            raise ValueError(f"Unknown dtype {dtype}")

        # ensure old data is gone
        if self.table_exists(table):
            self.conn.sql(f"""
                DELETE from {table}
                 WHERE name = '{name}'
                   AND exp_id = '{exp_id}'
            """)

        self.create_or_append(table=table, local_df=d)


    def uac_experiment_md(self,
                          expdict: dict) -> None:

        import pandas as pd

        expdata = pd.Series(expdict)
        ensure_str_columns = [
            'abstract', 'author', 'doi', 'dataset', 'experiment',
            'genes', 'organism', 'pubmed', 'study',
            'title', 'version', 'layer_type', 'layer_name']
        ensure_int_columns = ['year', 'dataset_id', 'full_experiment_id',
                              'study_id']

        for col in ensure_str_columns:
            if col not in expdata:
                expdata[col] = ''
            else:
                expdata[col] = str(expdata[col])

        for col in ensure_int_columns:
            if col not in expdata:
                expdata[col] = 0
            else:
                expdata[col] = int(expdata[col])

        return self.uac(
            table = 'experiment_md',
            local_df = expdata,
            ukey = 'dataset')


    def uac(self,
            table: str,
            local_df: Union["pd.DataFrame", "pd.Series", dict],
            ukey: str,
            ) -> None:
        """Check if the record exists - if so - update, otherwise create_or_append

        Name is from Update Append Create

        """

        import pandas as pd

        local_df = self._check_incoming_df(local_df)
        if not self.table_exists(table):
            return self.create_or_append(table, local_df)

        for _, row in local_df.iterrows():
            # delete (if exists)
            sql = f'''
                DELETE FROM {table}
                 WHERE "{ukey}" = '{row[ukey]}'
            '''
            self.conn.sql(sql)

        # & insert
        sql = f'''
            INSERT INTO '{table}'
            SELECT * FROM local_df
        '''

        self.conn.sql(sql)


    def create_or_append(
            self,
            table: str,
            local_df: Union["pd.DataFrame", "pd.Series", dict],
    ) -> None:

        import pandas as pd
        local_df = self._check_incoming_df(local_df)

        lg.debug(f"appending to {table} dataframe { local_df.shape }")
        local_df = local_df.sort_index(axis=1)

        #create a table - or if it exists - append
        if not self.table_exists(table):
            #lg("create & insert", table, local_df.shape)
            sql = f"CREATE TABLE '{table}' AS SELECT * FROM local_df"
            lg.debug(sql)
            self.sql(sql)
        else:
            #lg("append tbl", table, local_df.shape)
            sql = f"INSERT INTO '{table}' SELECT * FROM local_df"
            lg.debug(sql)
            self.sql(sql)
