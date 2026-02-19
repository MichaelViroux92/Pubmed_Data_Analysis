import streamlit as st
import requests 
import pandas as pd
import matplotlib.pyplot as plt

st.title("PubMed Search")

# Query search
search_term = st.text_input("Enter a search term:")

if st.button("Search"):
    response = requests.post(
        "http://fastapi-container:8000/search_query",
        json={"search_term": search_term}
    )
    if response.status_code == 200:
        st.json(response.json())
    else:
        st.error(f"Error: {response.status_code}")


# Trigger topic modeling
if st.button("Show subtopics"):
    response = requests.post(
        "http://fastapi-container:8000/clustering/topics/input"
    )

    if response.status_code == 200:
        data = response.json()
        topics = data.get("topics", [])

        st.write(f"Number of topics: {data.get('num_topics')}")

        for topic in topics:
            st.write(
                f"Topic {topic['Topic']} "
                f"({topic['Count']} docs): "
                f"{topic['Name']}"
            )
    else:
        st.error(f"Failed: {response.status_code}")


# Data analysis

# Generate a graph showing the number of articles per subtopic
if st.button("Number of articles per subtopic"):
    response = requests.get(  
        "http://fastapi-container:8000/labeled_data", 
    )

    if response.status_code == 200:
        data = response.json()
        df_data = pd.DataFrame(data["labeled_data"])

        subtopic_counts = df_data['topic_label'].value_counts()
        
        plt.figure(figsize=(10, 6))
        subtopic_counts.plot(kind='bar', color='skyblue')
        plt.title('Number of Articles per subtopic')
        plt.xlabel('Subtopic')
        plt.ylabel('Number of Articles')
        plt.xticks(rotation=45, ha='right')

        st.pyplot(plt)
    else:
        st.error(f"Failed  {response.status_code}")

