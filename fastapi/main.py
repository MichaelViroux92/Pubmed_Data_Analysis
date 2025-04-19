from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import JSONResponse

app = FastAPI()


# Searchquery endpoints    
@app.get("/search_query")
def search_query(search_term: str):
    response = requests.post(
        "http://working-container:5000/fetch_data",  # Replace with Docker service name
        json={"search_term": search_term}
    )
    return JSONResponse(content=response.json())



# Get inertia values for elbow curve
@app.get("/clustering/elbow/input")
def elbow_request():
    response = requests.post(
        "http://working-container:5000/clustering/elbow",
        json={}
    )
    return JSONResponse(content=response.json())


# Get the clustered labels
@app.get("/clustering/subtopics/input")
def get_clustering_subtopic_inputs(k: int = 20, num_words: int = 5):
    response = requests.post(
        "http://working-container:5000/clustering/subtopics",
        json={"k": k, "num_words": num_words}
    )
    return JSONResponse(content=response.json())




