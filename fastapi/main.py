from fastapi import FastAPI
from clustering import Clustering
from pydantic import BaseModel

app = FastAPI()
clustering = Clustering()
df = None # Global variable to store df

@app.get("/")
def read_root():
    return {"message": "API is working!"}


@app.post("/clustering/run")
def run_clustering():
    df = postgres.fetch_pubmed_data(["id", "abstract"])
    vectorized = clustering.vectorize(df)
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
    if df is None:
        df = postgres.fetch_pubmed_data(["id", "abstract"])
    vectorized = clustering.vectorize(df)
    clustering.fit_kmeans(vectorized, k=req.k)
    top_words = clustering.top5_cluster(num_words=req.num_words)
    return {"Top words": top_words}


#MVP
#Searchterm
# Streamlit -> Fastapi: input searchterm
# Fastapi - > pubmedapi and postgres functions
# Fastapi -> streamlit: return message, data loaded in database

#Elbowcurve
#Streamlit -> Fastapi: elbowcurve request
# Fastapi -> postgres load data and point to cluster elbow def
# Fastap -> streamlit: return inertia dict

#Clustering

# Streamlit -> Fastapi: input get clusterlabels. Default 20, user can choose number from elbow curve if generated. User can choose number of words. OK
# Fastapi -> postgres: load data and point to cluster defs. OK
# OPTIONAL. Extra api call to open AI for better cluster labels (function of clustering class)
# Fastapi -> sstreamlit: return cluster labels in dictionary. OK

# This is the Pydantic model for input validation


@app.get("/clustering/inputs")
def get_clustering_inputs(k: int = 20, num_words: int = 5):
    # Passing the default parameters to the POST endpoint
    return run_clustering(k=k, num_words=num_words)

class SubtopicRequest(BaseModel):
    k: int = 20  # Default value for k
    num_words: int = 5  # Default value for num_words

@app.post("/clustering/run")
def run_clustering(k: int, num_words: int):
    # Fetch the data, apply clustering, etc.
    df = postgres.fetch_pubmed_data(["id", "abstract"])
    vectorized = clustering.vectorize(df)
    inertia = clustering.elbowmethod(vectorized, kmax=30)
    return {"inertia": inertia}
