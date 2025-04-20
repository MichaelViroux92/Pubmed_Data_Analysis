from fastapi import FastAPI
from pydantic import BaseModel
from get_data.build_query import BuildQuery
from get_data.api_connection.pubmed_data_api import PubMedAPI
from db_connection.postgres import PostgresDB
from clustering.clustering import Clustering
import os
from dotenv import load_dotenv

load_dotenv()

app = FastAPI()
clustering = Clustering()
pubmedapi = PubMedAPI()
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
    full_query = BuildQuery.build_query(req.search_term)
    df = pubmedapi.fetch_data(full_query) 
    table = postgres.create_table()
    db = postgres.insert_pubmed_data(df)
    return {"message": "Data loaded into database"}


# Get inertia values for elbow curve
@app.post("/clustering/elbow")
def inertia_values():
    df_cluster = postgres.fetch_pubmed_data(["pmid", "abstract"])
    vectorized = clustering.vectorize(df_cluster)
    inertia = clustering.elbowmethod(vectorized, kmax=30)
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
    top_words = clustering.top5_cluster(num_words=req.num_words)
    return {"Top words": top_words}

