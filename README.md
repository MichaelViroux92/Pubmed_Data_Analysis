# PubMed Data Analysis Platform

A modular pipeline for fetching, clustering, and labeling scientific abstracts from PubMed using machine learning and OpenAI. Then this labeled data can be further analyzed for explorative purposes. The project uses a Docker-based microservice architecture and presents results in a Streamlit dashboard.

---

## Architecture Overview

This project uses **Docker Compose** to orchestrate the following services:

- **working-container**: 
  - Fetches PubMed data
  - Handles data processing & clustering logic
  - Interacts with PostgreSQL
- **fastapi-container**: 
  - API gateway
  - Exposes endpoints to trigger clustering and return results
- **streamlit-container**: 
  - Visualizes clustered data and insights
- **postgres-container**: 
  - Stores raw and processed data

---

## Features

- Fetch PubMed data and abstracts
- TF-IDF vectorization + KMeans clustering on abstracts
- Generate human-readable cluster labels using GPT (OpenAI API)
- Visualize labeled data in an interactive Streamlit dashboard
- Save & query results from PostgreSQL

---

## Getting Started

### 1. Clone the Repository

git clone https://github.com/MichaelViroux92/Pubmed_Data_Analysis.git

### 2. Environmental parameters

DB_USER=your_user
DB_PASSWORD=your_password
DB_NAME=your_database
DB_HOST=postgres-container
EMAIL=your_email
API_KEY=your_pubmed_api_key
OPENAI_API_KEY=your_open_ai_key


### 3. Build and start the application

docker-compose up --build

### 4. Access the Application

- Streamlit Dashboard: http://localhost:8501
- FastAPI UI: http://localhost:8001/docs

---

## Folder structure

```
Pubmed_Data_Analysis
│
├── data
│
├── fastapi
│   ├── Dockerfile
│   ├── main.py
│   └── requirements.txt
│
├── logic
│   ├── Dockerfile
│   ├── requirements.txt
│   ├── clustering
│   │   └── clustering.py
│   │
│   ├── db_connection
│   │   ├── __init__.py
│   │   └── postgres.py
│   │
│   └── get_data
│       ├── __pycache__/
│       ├── api_connection
│       │   ├── __pycache__/
│       │   ├── __init__.py
│       │   └── pubmed_data_api.py
│       │
│       ├── __init__.py
│       ├── build_query.py
│       └── api_endpoints.py
│
├── streamlit
│   ├── app.py
│   ├── Dockerfile
│   └── requirements.txt
│
├── .env
├── .gitignore
├── docker-compose.yml
└── README.md
```

---

## API Endpoints

### GET /search_query

- **Description**:
  - Triggers data fetching from PubMed using the provided search term. This sends a POST request to the working container.

- **Query Parameters**:
  - search_term (string): The search term used to query PubMed.

- **Returns**:
  - A success message once the data is fetched and loaded into the database.

### GET /clustering/elbow/input

- **Description**:
  - Fetches the inertia values used for plotting the elbow curve to determine the optimal number of clusters.

- **Returns**:
{
  "inertia": [value1, value2, ..., valueN]
}

### GET /clustering/subtopics/input

- **Description**:
  - Triggers clustering and labeling of articles into subtopics using K-Means.

- **Query Parameters**:
  - k (int): Number of clusters (default: 20)
  - num_words (int): Number of descriptive words per cluster (default: 5)

- **Returns**:
{
  "Labels": {
    "Cluster 0": "Vitamin D Metabolism",
    "Cluster 1": "Immune Function",
    ...
  }
}

### GET /labeled_data

- **Description**:
  - Fetches the labeled PubMed data from the Postgres database.

- **Returns**:
{
  "df_labeled_key": [
    {
      "pmid": "12345678",
      "title": "...",
      "abstract": "...",
      "cluster_name": "Immunology and Infection"
    },
    ...
  ]
}

### POST /fetch_data

- **Description**:
  - Fetches data from PubMed, stores it in the PostgreSQL database.

- **Request Body**:
{
  "search_term": "vitamin d"
}

- **Returns**:
{
  "message": "Data loaded into database"
}

### POST /clustering/elbow

- **Description**:
  - Computes inertia values from vectorized abstracts for elbow curve visualization.

- **Returns**:
{
  "inertia": [10.1, 8.4, 6.5, ...]
}

### POST /clustering/subtopics

- **Description**:
  - Performs K-Means clustering on article abstracts.
  - Uses openAI API to assign descriptive cluster names.

- **Request Body**:
{
  "k": 4,
  "num_words": 5
}

- **Returns**:
{
  "Labels": {
    "Cluster 0": "Bone Health",
    "Cluster 1": "Autoimmune Disease",
    ...
  }
}

### POST /labeled_dataframe

- **Description**:
  - Fetches the labeled DataFrame (as JSON records) from the database for further processing or visualization.

- **Returns**:
{
  "df_labeled_key": [
    {
      "pmid": "98765432",
      "title": "...",
      "abstract": "...",
      "cluster_name": "Metabolic Disorders"
    },
    ...
  ]
}

