import streamlit as st
import requests 
import pandas as pd
import matplotlib.pyplot as plt

st.title("PubMed Search")

# Query search
search_term = st.text_input("Enter a search term:")

if st.button("Search"):
    response = requests.get(
        "http://fastapi-container:8000/search_query",
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

        plt.figure(figsize=(8, 6))
        plt.plot(range(1, len(inertia_values) + 1), inertia_values, marker='o', linestyle='-', color='b')
        plt.title("Elbow Curve for KMeans Clustering")
        plt.xlabel("Number of Clusters (k)")
        plt.ylabel("Inertia")
        plt.grid(True)
    else:
        st.error(f"Failed to compute elbow curve: {response.status_code}")

# Clustered labels
# Streamlit input fields for k and num_words
k = st.slider("Number of clusters (k)", min_value=2, max_value=30, value=20)
num_words = st.slider("Number of top words per cluster", min_value=1, max_value=10, value=5)

# Streamlit button to trigger the clustering subtopics computation
if st.button("Show subtopics"):
    response = requests.get(  
        "http://fastapi-container:8000/clustering/subtopics/input", 
        params={"k": k, "num_words": num_words} 
    )

    if response.status_code == 200:
        labels = response.json().get('Labels', [])
        
        st.write("Labels for each cluster:")
        for i, label in labels.items():
            st.write(f"Subtopic {i}: {label}")
    else:
        st.error(f"Failed to retrieve clustering subtopics: {response.status_code}")


# Generate a graph showing the number of articles per subtopic
if st.button("Number of articles per subtopic"):
    response = requests.get(  
        "http://fastapi-container:8000/labeled_data", 
    )

    if response.status_code == 200:
        labeled_data = response.json()
        df_labeled = pd.DataFrame(labeled_data['df_labeled_key'])

        subtopic_counts = df_labeled['cluster_name'].value_counts()
        
        plt.figure(figsize=(10, 6))
        subtopic_counts.plot(kind='bar', color='skyblue')
        plt.title('Number of Articles per subtopic')
        plt.xlabel('Subtopic')
        plt.ylabel('Number of Articles')
        plt.xticks(rotation=45, ha='right')

        st.pyplot(plt)
    else:
        st.error(f"Failed  {response.status_code}")



