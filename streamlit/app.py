import streamlit as st


st.title("PubMed Search")

search_term = st.text_input("Enter a search term:")

if st.button("Search"):
    response = requests.get(
        "http://fastapi-container:8000/search_query",
        params={"search_term": search_term}
    )
    st.json(response.json())