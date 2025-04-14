class BuildQuery:
    @staticmethod
    def build_query(searchterm):
        article_type = ('[TIAB] AND ("Clinical Trial"[PT] OR "Randomized Controlled Trial"[PT] OR "Meta-Analysis"[PT] '
                        'OR "Systematic Review"[PT] OR "Comparative Study"[PT] OR "Observational Study"[PT])')

        full_query = searchterm + article_type
        print(f"Full query: {full_query}")
        return full_query
    
    