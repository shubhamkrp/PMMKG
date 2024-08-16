import requests
from bs4 import BeautifulSoup
import pandas as pd
from similarity_bw_pubmed_do import similar


# Define the base URL for the PubMed API
PUBMED_API_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Define the function to fetch disease terms related to the given symptom
def fetch_terms(symptom):
    # Replace spaces with '+' for URL encoding
    symptom_query = symptom.replace(" ", "+")

    # Construct the search query URL
    search_url = f"{PUBMED_API_BASE_URL}esearch.fcgi?db=pubmed&term={symptom_query}[MeSH Terms]&retmax=159"

    # Send the request to PubMed API
    response = requests.get(search_url)
    response.raise_for_status()

    soup = BeautifulSoup(response.content, 'xml')

    # Extract the PubMed IDs (PMIDs) from the response
    pmids = [id_tag.text for id_tag in soup.find_all('Id')]

    # Initialize a set to store fetched terms
    disease_terms = set()

    # Fetch details for each PMID
    for pmid in pmids:
        # Construct the fetch details URL
        fetch_url = f"{PUBMED_API_BASE_URL}efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"

        # Send the request to PubMed API
        fetch_response = requests.get(fetch_url)
        fetch_response.raise_for_status()

        # Parse the response XML using 'html.parser'
        fetch_soup = BeautifulSoup(fetch_response.content, 'xml')

        # Extract MeSH terms related to diseases
        mesh_headings = fetch_soup.find_all('MeshHeading')
        for mesh_heading in mesh_headings:
            descriptor = mesh_heading.DescriptorName
            if descriptor and descriptor.get('MajorTopicYN') == 'Y':
                disease_terms.add(descriptor.text)

    return disease_terms


symptom_file="mesh_symptoms.csv"
symptom_df=pd.read_csv(symptom_file, delimiter=';')

#select name column
column_name="Name"
symptom_term=symptom_df[column_name]

# print(symptom_term.head(2))

COLUMN_NAMES=["Symptom", "Related terms"]
df=pd.DataFrame(columns=COLUMN_NAMES)

no_df=pd.DataFrame(columns=["Symptom with no terms"])

for symptom in symptom_term:
    terms = fetch_terms(symptom)
    print(f"\nTerms associated with '{symptom}': {len(terms)}")
    if(len(terms)==0):
        no_df.loc[len(no_df)]=[symptom]
        continue
    term_string = ';'.join(terms)
    df.loc[len(df)]=[symptom,term_string]

no_df.to_csv("OUTPUT/mesh_symptom_with_no_terms.csv", sep=",", mode='a', index=False, header=False)
df.to_csv("OUTPUT/fetched_terms_for_mesh_symptoms.csv", sep=",", mode='a', index=False, header=False)
