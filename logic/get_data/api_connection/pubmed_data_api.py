from Bio import Entrez
import time
import json
import pandas as pd
from datetime import datetime
import os
from dotenv import load_dotenv
from get_data.build_query import BuildQuery

load_dotenv()

class PubMedAPI:
    def __init__(self):
        self.email = os.getenv("EMAIL")
        self.api_key = os.getenv("API_KEY")
        Entrez.email = self.email
        Entrez.api_key = self.api_key

    def fetch_data(self, full_query, count=9999):
        df = pd.DataFrame(columns=['PMID', 'Title', 'Abstract', 'Authors', 'Journal', 'Keywords', 'URL', 'Affiliations'])

        print(f"Email: {self.email}")
        print(f"API key: {self.api_key}")
        print(f"Query: {full_query}")

        # Fetching article IDs using esearch
        handle = Entrez.esearch(db="pubmed", term=full_query, retmax=min(count, 9999))
        search_results = Entrez.read(handle)
        pmids = search_results['IdList']

        # Fetching the actual article data using efetch
        handle = Entrez.efetch(db="pubmed", id=pmids, retmode="xml", api_key=Entrez.api_key)
        data = Entrez.read(handle)

        # Iterate over the fetched articles and extract details
        for record in data.get('PubmedArticle', []):
            pmid = record['MedlineCitation']['PMID']
            title = record['MedlineCitation']['Article']['ArticleTitle']

            abstract = ' '.join(record['MedlineCitation']['Article']['Abstract']['AbstractText']) if 'Abstract' in \
                                record['MedlineCitation']['Article'] and 'AbstractText' in \
                                record['MedlineCitation']['Article']['Abstract'] else ''

            authors = ', '.join(
                f"{author.get('LastName', '')} {author.get('ForeName', '')}".strip()
                for author in record['MedlineCitation']['Article'].get('AuthorList', [])
            )

            affiliations = []
            for author in record['MedlineCitation']['Article'].get('AuthorList', []):
                if 'AffiliationInfo' in author and author['AffiliationInfo']:
                    affiliations.append(author['AffiliationInfo'][0]['Affiliation'])
            affiliations = '; '.join(set(affiliations))

            journal = record['MedlineCitation']['Article']['Journal']['Title']
            keywords = ', '.join(
                keyword['DescriptorName'] for keyword in record['MedlineCitation'].get('MeshHeadingList', [])
            )
            url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"

            # Append the article data to the DataFrame
            df.loc[len(df)] = [pmid, title, abstract, authors, journal, keywords, url, affiliations]

        print("Data fetching complete.")
        return df

# Get search term and build query
searchterm = input("Enter query search term: ")
full_query = BuildQuery.build_query(searchterm)

# Create PubMedAPI instance and fetch data
pubmedapi = PubMedAPI()
df = pubmedapi.fetch_data(full_query)

# Display the first few rows
print(df.head())

# Clean up the search term for file naming
searchterm_cleaned = searchterm.replace('"', '')
searchterm_for_filename = searchterm_cleaned.replace(" ", "_")

# Define file paths
folder_path = os.path.join("..", "data")
os.makedirs(folder_path, exist_ok=True)
file_path = os.path.join(folder_path, f"{searchterm_for_filename}_results_api_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")

# Save results to a CSV file
df.to_csv(file_path, index=False)
print(f"Data saved to {file_path}")

# TODO: Explore how to fetch more than 10000 results using the API if possible
