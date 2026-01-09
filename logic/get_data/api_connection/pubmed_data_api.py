from Bio import Entrez
import time
import json
import pandas as pd
from datetime import datetime
import os
from dotenv import load_dotenv
from get_data.build_query import BuildQuery
from dataclasses import dataclass
from typing import List, Optional

load_dotenv()


@dataclass
class Article:
    """Class for keeping track of article attributes."""
    pmid: str
    title: str
    abstract: str
    publicationtype: List[str]
    authors: List[str]
    affiliations: List[str]
    journal: Optional[str]
    keywords: List[str]
    url: str
   

def parse_pubmed_record(record) -> Article:
    article = record['MedlineCitation']['Article']

    #pmid
    pmid = record['MedlineCitation']['PMID']
    
    #title
    title = article.get('ArticleTitle', '')
    
    #abstract
    abstract_parts = (article.get("Abstract", {}).get("AbstractText", []))
    processed_parts = []
    for part in abstract_parts:
        if isinstance(part, dict):
            label = part.get("@Label")
            text = part.get("#text", "")
            processed_parts.append(f"{label}: {text}" if label else text)
        else:
            processed_parts.append(str(part))
    abstract = " ".join(processed_parts)

    #publicationtype
    publicationtype = [
        pt.get('#text', '').strip() for pt in article.get('PublicationTypeList', [])
    ]

    #authors
    authors = [
        f"{author.get('LastName'), ''} {author.get('ForeName'), ''}".strip()
        for author in article.get('AuthorList', [])
        if author.get('LastName') or author.get('Forename')
    ]

    #affiliations
    affiliations = set()
    for author in article.get('AuthorList', []):
        if 'Affiliation' in author:
            affiliations.add(author['Affiliation'].strip())
        for aff_info in author.get('AffiliationInfo', []):
            if 'Affiliation' in aff_info:
                affiliations.add(aff_info['Affiliation'].strip())

    affiliations = list(affiliations)


    #journal
    journal = article.get('Journal', {}).get('Title', "")

    #keywords
    keywords = [
        keyword['DescriptorName'].get('#text', '') for keyword in record['MedlineCitation'].get('MeshHeadingList', [])
    ]
    
    #url
    url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}"

    return Article(
        pmid = pmid,
        title = title,
        abstract = abstract,
        publicationtype = publicationtype,
        authors = authors,
        affiliations = affiliations,
        journal = journal,
        keywords = keywords,
        url = url,
    )


def fetch_pubmed_articles(query: str, max_results: int = 9999) -> List[Article]:
    Entrez.email = os.getenv("EMAIL")
    Entrez.api_key = os.getenv("API_KEY")

    # Fetching article IDs using esearch
    handle = Entrez.esearch(db="pubmed", term=query, retmax=min(max_results, 9999))
    search_results = Entrez.read(handle)
    pmids = search_results['IdList']

    if not pmids:
        return []
    
    # Fetching the actual article data using efetch
    handle = Entrez.efetch(db="pubmed", id=pmids, retmode="xml", api_key=Entrez.api_key)
    data = Entrez.read(handle)

    articles = [
        parse_pubmed_record(record) for record in data.get('PubmedArticle', [])
    ]

    return articles



'''

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

'''