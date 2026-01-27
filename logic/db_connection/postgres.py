import psycopg2
from psycopg2.extras import execute_values
import pandas as pd
from typing import List, Optional
from get_data.api_connection.pubmed_data_api import Article

class PostgresDB:
    def __init__(self, host: str, dbname: str, user: str, password: str, port: int= 5432):
        try:    
            self.conn = psycopg2.connect(
                host=host,
                dbname=dbname,
                user=user,
                password=password,
                port=port
            )
        except psycopg2.Error as e:
            raise ConnectionError(f"Error connecting to PostgreSQL: {e}")
        
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def create_table(self) -> None:
        query = """
            CREATE TABLE IF NOT EXISTS pubmed_articles(
                pmid TEXT PRIMARY KEY,
                title TEXT,
                abstract TEXT,
                publicationtype TEXT,
                authors TEXT,
                affiliations TEXT,
                journal TEXT,
                keywords TEXT,
                url TEXT,
                cluster_name TEXT
            );
            """
        try:
            with self.conn.cursor() as cursor:
                cursor.execute(query)
            self.conn.commit()
        except psycopg2.Error as e:
            self.conn.rollback()
            raise RuntimeError(f"Failed to create table: {e}")

    def insert_articles(self, articles: List[Article]) -> None:
        if not articles:
            return
        query = """
        INSERT INTO pubmed_articles
        (pmid, title, abstract, publicationtype, authors, affiliations, journal, keywords, url)
        VALUES %s
        ON CONFLICT (pmid) DO NOTHING;
        """
        values = [
            (
                a.pmid,
                a.title,
                a.abstract,
                ", ".join(a.publicationtype),
                ", ".join(a.authors),
                "; ".join(a.affiliations),
                a.journal,
                ", ".join(a.keywords),
                a.url,
            )
            for a in articles
        ]
     
        try:
            with self.conn.cursor() as cursor:
                execute_values(cursor, query, values)
            self.conn.commit()
        except psycopg2.Error as e:
            self.conn.rollback()
            raise RuntimeError(f"Failed to insert articles: {e}")

    def fetch_articles_df(self, columns: Optional[List[str]] = None) -> pd.DataFrame:
        if columns is None:
            query = "SELECT * FROM pubmed_articles"
        else: 
            query = f"SELECT {', '.join(columns)} FROM pubmed_articles"

        try:
            with self.conn.cursor() as cursor: 
                cursor.execute(query)
                rows = cursor.fetchall()
        
                if columns is None:
                    columns = [desc[0] for desc in cursor.description] # cursor.description: [('id',), ('name',)] tuple. [0] is column name

        except psycopg2.Error as e:
            raise RuntimeError(f"Failed to fetch data: {e}")

        # Build dataframe        
        return pd.DataFrame(rows, columns=columns)

    def update_cluster_labels(self, updates: List[tuple[str, str]]) -> None:
        if not updates:
            return
        
        query = """
        UPDATE pubmed_articles as p
        SET cluster_name = data.cluster_name
        FROM (VALUES %s) AS data(pmid, cluster_name)
        WHERE p.pmid = data.pmid
        """
        try:
            with self.conn.cursor() as cursor:
                execute_values(cursor, query, updates)
            self.conn.commit()
        except psycopg2.Error as e:
            self.conn.rollback()
            raise RuntimeError(f"Failed to bulk update cluster labels: {e}")

    def close(self) -> None:
        if self.conn:
            self.conn.close()
