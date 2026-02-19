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
                topic_label TEXT,
                topic_id INTEGER
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

    def _fetch_articles(self, columns: Optional[List[str]] = None) -> Tuple[List[str], List[tuple]]:
        if columns is None:
            query = "SELECT * FROM pubmed_articles"
        else: 
            query = f"SELECT {', '.join(columns)} FROM pubmed_articles"

        try:
            with self.conn.cursor() as cursor: 
                cursor.execute(query)
                rows = cursor.fetchall()
                colnames = [desc[0] for desc in cursor.description]

        except psycopg2.Error as e:
            raise RuntimeError(f"Failed to fetch data: {e}")
              
        return colnames, rows
    
    def fetch_articles_ML(self, columns: List[str]) -> Tuple[List, ...]:
        _, rows = self._fetch_articles(columns)

        if not rows:
            return tuple([] for _ in columns)
        
        return tuple(map(list, zip(*rows))) #for each element inside zip(*rows), apply list() to it

    
    def fetch_articles_JSON(self, columns: List[str]) -> List[dict]:
        colnames, rows = self._fetch_articles(columns)

        return [dict(zip(colnames, row)) for row in rows]

    def update_cluster_labels(self, updates: List[tuple[str, str, int]]) -> None:
        if not updates:
            return
        
        query = """
        UPDATE pubmed_articles as p
        SET 
            topic_label = data.topic_label
            topic_id = data.topic_id
        FROM (VALUES %s) AS data(pmid, topic_label, topic_id)
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
