
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans


class Clustering:
    def __init__(self):
        self.vectorizer_abstracts = TfidfVectorizer(stop_words='english', ngram_range=(1,2))
        self.inertia = []
        self.kmeans_model = None
        self.vectorized_abstracts = None
        self.df = None 

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

    def top5_cluster(self, num_words=5):
        feature_names_abstracts = self.vectorizer_abstracts.get_feature_names_out()
        top_words_per_cluster = {}
        
        for i, centroid in enumerate(self.kmeans_model.cluster_centers_):
            top_indices = centroid.argsort()[-num_words:][::-1]
            top_words = [feature_names_abstracts[idx] for idx in top_indices]
            top_words_per_cluster[i] = top_words

        return top_words_per_cluster

