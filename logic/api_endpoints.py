from fastapi import FastAPI, Depends
from pydantic import BaseModel
from get_data.build_query import build_query
from get_data.api_connection.pubmed_data_api import fetch_pubmed_articles
from db_connection.postgres import PostgresDB
import os
from dotenv import load_dotenv
import clustering.topic_modeling as tm
from clustering.config import EMBEDDING_MODEL

load_dotenv()

app = FastAPI()

def get_db():
    db = PostgresDB(
        host=os.getenv("DB_HOST"),
        dbname=os.getenv("DB_NAME"),
        user=os.getenv("DB_USER"),
        password=os.getenv("DB_PASSWORD")
    )
    try:
        yield db
    finally:
        db.close()


# Searchterm query to pubmed
class FetchData(BaseModel):
    search_term: str

@app.post("/fetch_data")
def fetch_data(req: FetchData, db: PostgresDB = Depends(get_db)):
    full_query = build_query(req.search_term)
    articles = fetch_pubmed_articles(full_query) 
    db.create_table()
    db.insert_articles(articles)
    return {"message": "Data loaded into database"}

# Topic labeling
@app.post("/clustering/topics")
def topic_labeling(db: PostgresDB = Depends(get_db)):
    pmids, abstracts = db.fetch_articles_ML(["pmid", "abstract"])
    embedding_model = tm.build_embedder(EMBEDDING_MODEL)
    embedded_abstracts = tm.encode_abstracts(embedding_model, abstracts)
    topic_model = tm.build_topic_model()
    topics, _ = topic_model.fit_transform(abstracts, embedded_abstracts)
    representation_model = tm.build_represenation_model()
    topic_model.update_topics(abstracts, representation_model=representation_model)
    topic_info = topic_model.get_topic_info()
    topic_id_to_label = dict(zip(topic_info.Topic, topic_info.Name))
    update_db_labels = tm.map_pmid_to_topic(
        pmids=pmids,
        topics=topics,
        topic_id_to_label=topic_id_to_label,
    )
    db.update_cluster_labels(update_db_labels)
    return {
        "num_topics": len(topic_info) - 1,
        "topics": topic_info.to_dict(orient='records')
    }


# Fetch data for streamlit analysis
@app.get("/labeled_data")
def get_labeled_dataframe(db: PostgresDB = Depends(get_db)):
    data = db.fetch_articles_JSON()
    return {"labeled_data": data}



# ----------------------------------------------------------------------

"""
# Get inertia values for elbow curve
@app.post("/clustering/elbow")
def inertia_values():
    df_cluster = postgres.fetch_pubmed_data(["pmid", "abstract"])
    vectorized = clustering.vectorize(df_cluster)
    inertia = clustering.elbowmethod(vectorized, kmax=25)
    return {"inertia": inertia}



# Get the clustered labels
class SubtopicRequest(BaseModel): # The class is meant for the post endpoint
    k: int = 20
    num_words: int = 5


@app.post("/clustering/subtopics")
def clustering_subtopics(req: SubtopicRequest):
    df_cluster = postgres.fetch_pubmed_data(["pmid", "abstract"])
    vectorized = clustering.vectorize(df_cluster)
    clustering.fit_kmeans(vectorized, k=req.k)

    descriptive_cluster_names = clustering.descriptive_names_clusters(num_words=req.num_words)
    pmid_to_name = clustering.pmid_clustername_mapping(descriptive_cluster_names)

    for pmid, cluster_name in pmid_to_name.items():
        postgres.update_cluster_labels(pmid, cluster_name)

    return {"Labels": descriptive_cluster_names}

"""

