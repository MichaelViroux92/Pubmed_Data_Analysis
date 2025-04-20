# PubMed Data Analysis Platform

A modular pipeline for fetching, clustering, and labeling scientific abstracts from PubMed using machine learning and OpenAI. Then this labeled data can be further analyzed for explorative purposes. The project uses a **Docker-based microservice architecture** and presents results in a **Streamlit dashboard**.

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

## 🚀 Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/pubmed-clustering.git
cd pubmed-clustering

### 2. Environmental parameters

OPENAI_API_KEY=your_openai_key
POSTGRES_USER=your_user
POSTGRES_PASSWORD=your_password
POSTGRES_DB=your_database


### 3. Build and start the application

docker-compose up --build

