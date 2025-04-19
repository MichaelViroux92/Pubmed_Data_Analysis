from fastapi import FastAPI
from clustering import Clustering
from pydantic import BaseModel


app = FastAPI()
clustering = Clustering()
pubmedapi = PubMedAPI()
postgres = PostgresDB()

global df_cluster = None # Global variable to store df

@app.get("/")
def read_root():
    return {"message": "API is working!"}


# Searchquery endpoints
class FetchData(BaseModel):
    search_term: str

@app.get("/search_query/input")
def search_query(search_term: str):
    return fetch_data()
    

@app.post("/fetch_data")
def get_data(req: FetchData)
    full_query = BuildQuery.build_query(req.search_term)
    df = pubmed_data_api.fetch_data(full_query) 
    table = postgres.create_table()
    db = postgres.insert_pubmed_data(df)
    return {"message": "Data loaded into database"}


# Get inertia values for elbow curve
@app.get("/clustering/elbow/request")
def elbow_request():
    return inertia_values()

@app.post("/clustering/elbow")
def inertia_values():
    df_cluster = postgres.fetch_pubmed_data(["id", "abstract"])
    vectorized = clustering.vectorize(df_cluster)
    inertia = clustering.elbowmethod(vectorized, kmax=30)
    return {"inertia": inertia}


# Get the clustered labels
class SubtopicRequest(BaseModel): # The class is meant for the post endpoint
    k: int = 20
    num_words: int = 5

@app.get("/clustering/subtopics/inputs")
def get_clustering_subtopic_inputs(k: int = 20, num_words: int = 5):
    # Passing the default parameters to the POST endpoint
    return clustering_subtopics(k=k, num_words=num_words)


@app.post("/clustering/subtopics")
def clustering_subtopics(req: SubtopicRequest):
    if df_cluster is None:
        df_cluster = postgres.fetch_pubmed_data(["id", "abstract"])
    vectorized = clustering.vectorize(df_cluster)
    clustering.fit_kmeans(vectorized, k=req.k)
    top_words = clustering.top5_cluster(num_words=req.num_words)
    return {"Top words": top_words}


