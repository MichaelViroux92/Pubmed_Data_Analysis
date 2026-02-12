# Embedding model
EMBEDDING_MODEL = "pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb"

# UMAP parameters
UMAP_CONFIG = {
    n_neighbors: 30,
    n_components: 8,
    min_dist: 0.0,
    metric: 'cosine',
    random_state: 42
}

# HDBSCAN parameters
HDBSCAN_CONFIG = {
    min_cluster_size: 30,
    min_samples: 15,
    metric: "euclidean",
    cluster_selection_method: "eom",
}

# CountVectorizer parameters
VECTORIZER_CONFIG = {
    ngram_range: (1, 2),
    min_df: 5,
    max_df: 0.95,
    stop_words: "english"
}

# OpenAI
OPENAI_MODEL = "gpt-4o-mini"