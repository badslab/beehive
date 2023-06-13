

import streamlit as st
from termite import db


def welcome():
    st.title("Welcome!")
    exp_data = db.raw_sql('select * from help_experiments')
    st.write(exp_data)
