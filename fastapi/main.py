from fastapi import FastAPI
from clustering import Clustering
from pydantic import BaseModel

app = FastAPI()
clustering = Clustering()

@app.get("/")
def read_root():
    return {"message": "API is working!"}


@app.post("/clustering/run")
def run_clustering():
    df = postgres.fetch_pubmed_data(["id", "abstract"])
    vectorized = clustering.vectorize(df)
    inertia = clustering.elbowmethod(vectorized, kmax=30)
    return {"inertia": inertia}


class SubtopicRequest(BaseModel):
    k: int = 20
    num_words: int = 5

@app.post("/clustering/subtopics")
def subtopics(req: SubtopicRequest):
    df = postgres.fetch_pubmed_data(["id", "abstract"])
    vectorized = clustering.vectorize(df)
    clustering.fit_kmeans(vectorized, k=req.k)
    top_words = clustering.top5_cluster(num_words=req.num_words)