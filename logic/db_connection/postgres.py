import psycopg2
from psycopg2.extras import execute_values
import pandas as pd

class PostgresDB:
    def __init__(self, host, dbname, user, password, port=5432):
        self.conn = psycopg2.connect(
            host=host,
            dbname=dbname,
            user=user,
            password=password,
            port=port
        )
        self.cursor = self.conn.cursor()

# pmid, title, abstract, authors, journal, keywords, url, affiliations

    def create_table(self):
        self.cursor.execute("""
                            CREATE TABLE IF NOT EXISTS pubmed_articles
                            (
                            pmid TEXT PRIMARY KEY,
                            title TEXT,
                            abstract TEXT,
                            authors TEXT,
                            journal TEXT,
                            keywords TEXT,
                            url TEXT,
                            affiliations TEXT
                            );
                            """)
        self.conn.commit()

    def insert_pubmed_data(self, df):
        query = """
        INSERT INTO pubmed_articles
        (pmid, title, abstract, authors, journal, keywords, url, affiliations)
        VALUES %s
        ON CONFLICT (pmid) DO NOTHING;
        """
        values = list(df.itertuples(index=False, name=None))
        execute_values(self.cursor, query, values)
        self.conn.commit()

    def fetch_pubmed_data(self, columns=None):
        if columns is None:
            query = "SELECT * FROM pubmed_articles"
        else:
            selected_columns = ", ".join(columns)
            query = f"SELECT {selected_columns} FROM pubmed_articles"

        self.cursor.execute(query)
        rows = self.cursor.fetchall()

        # Build dataframe
        if columns is None:
            columns = [desc[0] for desc in self.cursor.description] # cursor.description: [('id',), ('name',)]. So we just get all columns

        df = pd.DataFrame(rows, columns=columns)
        return df

    def close(self):
        self.cursor.close()
        self.conn.close()