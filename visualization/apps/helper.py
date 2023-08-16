
import pandas as pd
import streamlit as st

from termite.vis import util
from termite.db import get_conn


def tablepeek():
    db = get_conn()
        
    rv = {}
    for _, table in db.sql('SHOW ALL TABLES').df().iterrows():
        t = table['name']
        c = db.sql(f"SELECT count(*) FROM {t}").df().iloc[0,0]
        rv[t] = c
    
    tablecount = pd.Series(rv)
    tablecount.columns = 'No Records'
    table = st.sidebar.selectbox("table", tablecount.index)
    head = util.execute_sql(f"select * from {table} limit 5")
    st.dataframe(tablecount, use_container_width=True)
    if st.sidebar.checkbox('transpose?'):
        head = head.T
    st.dataframe(head, use_container_width=True)


def raw_sql():
    sql = st.text_area(
        "Sql", height=10, 
        value = "SELECT * FROM expr LIMIT 5", key='sql')
    
    res = util.execute_sql(sql)
    st.write(res)

    
