import streamlit as st
import requests  # Make sure to import requests

st.title("PubMed Search")

# Query search
search_term = st.text_input("Enter a search term:")

if st.button("Search"):
    response = requests.get(
        "http://fastapi-container:8000/search_query",  # Correct URL for FastAPI container
        params={"search_term": search_term}
    )
    if response.status_code == 200:
        st.json(response.json())
    else:
        st.error(f"Error: {response.status_code}")

# Elbow curve
if st.button("Compute Elbow Curve"):
    response = requests.get("http://fastapi-container:8000/clustering/elbow/input")
    
    if response.status_code == 200:
        inertia_data = response.json()
        inertia_values = inertia_data['inertia']
        
        st.write("Inertia values:", inertia_values)
    else:
        st.error(f"Failed to compute elbow curve: {response.status_code}")

# Clustered labels
# Streamlit input fields for k and num_words
k = st.slider("Number of clusters (k)", min_value=2, max_value=30, value=20)
num_words = st.slider("Number of top words per cluster", min_value=1, max_value=10, value=5)

# Streamlit button to trigger the clustering subtopics computation
if st.button("Get Clustering Subtopics"):
    response = requests.get(  # Use GET here, as per your original request
        "http://fastapi-container:8000/clustering/subtopics/input",  # Correct URL
        params={"k": k, "num_words": num_words}  # Changed to use `params` for GET request
    )

    if response.status_code == 200:
        labels = response.json().get('Labels', [])
        
        # Display the top words for each cluster
        st.write("Labels for each cluster:")
        for i, label in labels.items():
            st.write(f"Subtopic {i}: {label}")
    else:
        st.error(f"Failed to retrieve clustering subtopics: {response.status_code}")

