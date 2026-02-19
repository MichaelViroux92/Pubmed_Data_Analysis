from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import JSONResponse
import requests

app = FastAPI()

# Searchquery endpoints

class FetchData(BaseModel):
    search_term: str

@app.post("/search_query")
def search_query(req: FetchData):
    response = requests.post(
        "http://working-container:8001/fetch_data",  # Replace with Docker service name
        json={"search_term": req.search_term}
    )
    return JSONResponse(content=response.json())


# Get the generated topics
@app.post("/clustering/topics/input")
def get_clustering_subtopic_inputs():
    response = requests.post(
        "http://working-container:8001/clustering/topics",
        json={}
    )
    return JSONResponse(content=response.json())


# Get data for analysis
@app.get("/labeled_data")
def get_labeled_data():
    response = requests.get(
        "http://working-container:8001/labeled_data"
    )
    return JSONResponse(content=response.json())




# ------------------------------------------------------------------------

"""
# Get inertia values for elbow curve
@app.get("/clustering/elbow/input")
def elbow_request():
    response = requests.post(
        "http://working-container:8001/clustering/elbow",
        json={}
    )
    return JSONResponse(content=response.json())

"""