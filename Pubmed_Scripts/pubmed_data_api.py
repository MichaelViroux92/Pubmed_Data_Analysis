from Bio import Entrez
import time
import json
import pandas as pd


def build_query():
    searchterm = input("Enter query search term: ")
    article_type = ('[IT] AND ("Clinical Trial"[PT] OR "Randomized Controlled Trial"[PT] OR "Meta-Analysis"[PT] '
                    'OR "Systematic Review"[PT] OR "Comparative Study"[PT] OR "Observational Study"[PT] '
                    'OR "Validation Study"[PT] OR "Case Reports"[PT] OR "Review"[PT])')

    full_query = searchterm + article_type
    return full_query


def history(full_query):
    API_key = input("Enter API_key:")

    # Get number of data records
    with Entrez.esearch(db="pubmed", term=query, retmax=1, api_key=API_key, usehistory="y") as handle:
        results = Entrez.read(handle)
        count = int(results["Count"])
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]

    print(f"Total articles found: {count}")
    print(f"WebEnv: {webenv}, QueryKey: {query_key}")

    return webenv, query_key, count, API_key


def fetch_data(webenv, query_key, count, API_key, batchsize=5000):
    max_records = min(count, 9999)
    df = pd.DataFrame(columns=['PMID', 'Title', 'Abstract', 'Authors', 'Journal', 'Keywords', 'URL', 'Affiliations'])

    for start in range(0, max_records, batchsize):
        print(f"Records {start + 1} to {min(start + batchsize, count)}...")

        # Request data from PubMed using efetch, and specify 'json' for retmode
        with Entrez.efetch(db="pubmed", query_key=query_key, webenv=webenv, retstart=start, retmax=batchsize,
                           retmode="xml", api_key=API_key) as handle:
            data = Entrez.read(handle)

            # Process each article
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

                # Append the data to the DataFrame
                df.loc[len(df)] = [pmid, title, abstract, authors, journal, keywords, url, affiliations]

        time.sleep(10)

    print("Data fetching complete.")
    return df

Entrez.email = input("Enter emailadres:")
query = build_query()
webenv, query_key, count, API_key = history(query)


if count > 0:
    df = fetch_data(webenv, query_key, count, API_key)
    print(df.head())  # Display a preview of the data
    df.to_csv("pubmed_results2.csv", index=False)  # Save results to a CSV file
    print("Data saved to pubmed_results.csv")
else:
    print("No articles found.")