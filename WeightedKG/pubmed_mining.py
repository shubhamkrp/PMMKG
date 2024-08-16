from Bio import Entrez
import csv, os

# Function to search PubMed
def search_pubmed(term, max_results=10):
    Entrez.email = "youemail@example.com"
    handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# Function to fetch article details
def fetch_details(id_list):
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

# Function to extract MeSH terms
def extract_mesh_terms(records):
    mesh_terms = set()
    for record in records["PubmedArticle"]:
        if "MeshHeadingList" in record["MedlineCitation"]:
            for mesh_heading in record["MedlineCitation"]["MeshHeadingList"]:
                mesh_terms.add(mesh_heading["DescriptorName"])
    return mesh_terms

# Function to extract title, abstract, and full text link
def extract_article_details(records):
    articles = []
    mesh_terms = set()
    for record in records["PubmedArticle"]:
        if "MeshHeadingList" in record["MedlineCitation"]:
            for mesh_heading in record["MedlineCitation"]["MeshHeadingList"]:
                mesh_terms.add(mesh_heading["DescriptorName"].lower())
        pubmed_id=record["MedlineCitation"]["PMID"]
        article = record["MedlineCitation"]["Article"]
        title = article.get("ArticleTitle", "No title available")
        abstract = article.get("Abstract", {}).get("AbstractText", ["No abstract available"])[0]
        full_text_link = None

        # Try to find full text link from PubmedData
        if "PubmedData" in record and "ArticleIdList" in record["PubmedData"]:
            for article_id in record["PubmedData"]["ArticleIdList"]:
                if article_id.attributes["IdType"] == "pmc":
                    full_text_link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{article_id}/"

        articles.append({
            "pubmed_id": pubmed_id,
            "title": title.lower(),
            "abstract": abstract.lower(),
            "full_text_link": full_text_link,
        })
    return mesh_terms, articles

def save_to_csv(articles, filename):
    keys=["pubmed_id", "title", "abstract", "full_text_link"]
    with open(filename, 'w', newline='', encoding='utf-8') as output_file: 
        dict_writer = csv.DictWriter(output_file, fieldnames=keys) 
        dict_writer.writeheader() 
        dict_writer.writerows(articles)


import pandas as pd
symp_file="mesh_symptoms.csv"
symp_df=pd.read_csv(symp_file, delimiter=';')
symp_term=symp_df["Name"]
print(symp_term.head())

# term = "Prosopagnosia"
for term in symp_term:
    term=term.lower()
    ids = search_pubmed(term, max_results=999999)
    details = fetch_details(ids)
    mesh_terms, articles = extract_article_details(details)
    print(len(mesh_terms))

    print(len(articles))
    term_name = "_".join(term.split(" "))

    if not os.path.exists(term_name):
        os.mkdir(term_name)

    with open(os.path.join(term_name, f'{term_name}.txt'), 'w') as f:
        f.write(f"{term}\n")
        for line in mesh_terms:
            f.write(f"{line.lower()}\n")

    save_to_csv(articles, os.path.join(term_name, f'{term_name}_articles.csv'))
