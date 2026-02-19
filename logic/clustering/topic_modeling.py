from clustering.config import (
    EMBEDDING_MODEL,
    UMAP_CONFIG,
    HDBSCAN_CONFIG,
    VECTORIZER_CONFIG,
    OPENAI_MODEL
)
from typing import List, Dict, Tuple
from sentence_transformers import SentenceTransformer
from umap import UMAP
from hdbscan import HDBSCAN
from bertopic import BERTopic
from sklearn.feature_extraction.text import CountVectorizer
from dotenv import load_dotenv
import os
import openai
from bertopic.representation import OpenAI
from clustering.prompts import TOPIC_LABEL_PROMPT

load_dotenv()
open_api_key = os.getenv("OPENAI_API_KEY")

# Embeddings
def build_embedder(model_name: str):
    return SentenceTransformer(model_name)

def encode_abstracts(model, abstracts: List[str]):
    return model.encode(abstracts)

# Modeling
def build_topic_model():
    return BERTopic(
        umap_model = UMAP(**UMAP_CONFIG),
        hdbscan_model = HDBSCAN(**HDBSCAN_CONFIG),
        vectorizer_model = CountVectorizer(**VECTORIZER_CONFIG)
    )
        
# Representation
def build_represenation_model():
    client = openai.OpenAI(api_key=open_api_key) #client, connections as seperate variable
    
    return OpenAI( #config values get passed inline
        client, 
        model=OPENAI_MODEL, 
        exponential_backoff=True, 
        prompt=TOPIC_LABEL_PROMPT
    )

# Mapping topics to articles
def map_pmid_to_topic(pmids: List[str], topics: List[int], topic_id_to_label: Dict[int, str]) -> List[Tuple[str, str]]:
    return [
         (pmid, topic_id_to_label.get(topic_id), topic_id)
         for pmid, topic_id in zip(pmids, topics)
    ]
