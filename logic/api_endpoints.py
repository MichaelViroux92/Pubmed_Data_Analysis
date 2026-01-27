from fastapi import FastAPI
from pydantic import BaseModel
from get_data.build_query import build_query
from get_data.api_connection.pubmed_data_api import fetch_pubmed_articles
from db_connection.postgres import PostgresDB
from clustering.clustering import Clustering
import os
from dotenv import load_dotenv

load_dotenv()

app = FastAPI()
clustering = Clustering()
postgres = PostgresDB(
    host=os.getenv("DB_HOST"),
    dbname=os.getenv("DB_NAME"),
    user=os.getenv("DB_USER"),
    password=os.getenv("DB_PASSWORD")
)


# Searchquery endpoints
class FetchData(BaseModel):
    search_term: str

@app.post("/fetch_data")
def fetch_data(req: FetchData):
    full_query = build_query(req.search_term)
    articles = fetch_pubmed_articles(full_query) 
    postgres.create_table()
    postgres.insert_articles(articles)
    return {"message": "Data loaded into database"}


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

# Get clustered dataframe
@app.post("/labeled_dataframe")
def get_labeled_dataframe():
    df_labeled = postgres.fetch_pubmed_data()
    df_labeled_json = df_labeled.to_dict(orient="records")
    return {"df_labeled_key": df_labeled_json}