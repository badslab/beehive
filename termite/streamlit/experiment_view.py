

import streamlit as st
from termite import db


def experiment_view():

    studies = db.raw_sql(
        '''SELECT study, max(year) as maxyear
             FROM experiment_md
         GROUP BY study
         ORDER BY maxyear''')
    
    study = st.selectbox("Study", 
        options=['Zoom in on a study'] +  list(studies['study']),
        label_visibility='collapsed', )
    search = st.text_input("Search string").strip()

    where = []
    if study != 'Zoom in on a study':
        where.append(f"study = '{study}'")
    if search:
        where.append(f"title LIKE '%{search}%'")

    sql = "SELECT * FROM experiment_md "
    
    if where:
        sql += " WHERE " + " AND ".join(where) 

    exp_data = db.raw_sql(sql)
    if len(exp_data) == 0:
        st.write("No experiments found")
        st.stop()
        
    # st.write(exp_data.head(2).T)

    exp_data['Title'] = exp_data.apply(
        lambda r: f"<b>{r['title']}</b>. <i>{r['author']}</i> ({r['year']})",
        axis=1)

    del exp_data['title']
    del exp_data['author']
    del exp_data['year']

    
    def create_links(row):
        # target _blank to open new window
        rv = ''
        if row['doi'] != 'unknown':
            rv += '<a target="_blank" href="https://doi.org/{}">D</a>'\
                .format(row['doi'].split()[0])
        if row['pubmed']:
            rv += " <a target='_blank' href=;https://pubmed.ncbi.nlm.nih.gov/{}'>P</a>"
        return rv
    
    exp_data['link'] = exp_data.apply(create_links, axis=1)
    
    exp_data = exp_data.reindex(
        ['Title', 'experiment', 'link'], axis=1)

    st.subheader("Experiments:")
    styled = exp_data.style
    st.markdown(styled.hide().to_html(), unsafe_allow_html=True)
    
