
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from dotenv import load_dotenv
import os
import openai
 
load_dotenv()

class Clustering:
    def __init__(self):
        self.vectorizer_abstracts = TfidfVectorizer(stop_words='english', ngram_range=(1,2))
        self.inertia = []
        self.kmeans_model = None
        self.vectorized_abstracts = None
        self.df = None
        self.openai_api_key = os.getenv("OPENAI_API_KEY")
        openai.api_key = self.openai_api_key  

    def vectorize(self, df):
        self.df = df
        df['abstract'] = df['abstract'].fillna("No abstract available")
        self.vectorized_abstracts = self.vectorizer_abstracts.fit_transform(df['abstract'])

        return self.vectorized_abstracts
    
    def elbowmethod(self, vectorized_abstracts, kmax):
        self.inertia.clear()

        k_values = range(2, kmax)
        for k in k_values:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            kmeans.fit(vectorized_abstracts)
            self.inertia.append(kmeans.inertia_)  # Store inertia (sum of squared distances)

        return self.inertia
    
    def fit_kmeans(self, vectorized_abstracts, k=20):
        self.kmeans_model = KMeans(n_clusters=k, random_state=42, n_init=10)
        self.kmeans_model.fit(vectorized_abstracts)

    def topwords_cluster(self, num_words):
        feature_names_abstracts = self.vectorizer_abstracts.get_feature_names_out()
        top_words_per_cluster = {}
        
        for i, centroid in enumerate(self.kmeans_model.cluster_centers_):
            top_indices = centroid.argsort()[-num_words:][::-1]
            top_words = [feature_names_abstracts[idx] for idx in top_indices]
            top_words_per_cluster[i] = top_words

        return top_words_per_cluster


    def get_cluster_name(self, top_words):
        prompt = f"Based on the following top words, generate a descriptive name for this cluster: {', '.join(top_words)}"
        
        response = openai.chat.completions.create(
            model="gpt-4o",
            messages=[
            {"role": "user", "content": prompt}
            ],
            max_tokens=10,
            n=1,
            stop=None,
            temperature=0.5
        )

        cluster_name = response.choices[0].message.content.strip()
        return cluster_name
    

    def descriptive_names_clusters(self, num_words=5):
        top_words_per_cluster = self.topwords_cluster(num_words=num_words)
        cluster_names = {}

        for cluster_id, top_words in top_words_per_cluster.items():
            cluster_name = self.get_cluster_name(top_words)
            cluster_names[cluster_id] = cluster_name

        return cluster_names