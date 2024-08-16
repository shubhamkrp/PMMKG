import requests
from bs4 import BeautifulSoup

# Define the base URL for the PubMed API
PUBMED_API_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Define the disease for which we want to fetch related symptom terms
disease = "fever"  # Replace with your desired disease

# Define the function to fetch symptom terms related to the given disease
def fetch_symptom_terms(disease):
    # Replace spaces with '+' for URL encoding
    disease_query = disease.replace(" ", "+")

    # Construct the search query URL
    search_url = f"{PUBMED_API_BASE_URL}esearch.fcgi?db=pubmed&term={disease_query}[MeSH Terms][Majr:NoExp]&retmax=10"

    # Send the request to PubMed API
    response = requests.get(search_url)
    response.raise_for_status()  # Raise an error if the request failed

    soup = BeautifulSoup(response.content, 'xml')

    # Extract the PubMed IDs (PMIDs) from the response
    pmids = [id_tag.text for id_tag in soup.find_all('Id')]

    # Initialize a set to store unique symptom terms
    symptom_terms = set()

    # Fetch details for each PMID
    for pmid in pmids:
        # Construct the fetch details URL
        fetch_url = f"{PUBMED_API_BASE_URL}efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"

        # Send the request to PubMed API
        fetch_response = requests.get(fetch_url)
        fetch_response.raise_for_status()

        # Parse the response XML using 'html.parser'
        fetch_soup = BeautifulSoup(fetch_response.content, 'xml')

        # Extract MeSH terms related to symptoms
        mesh_headings = fetch_soup.find_all('MeshHeading')
        for mesh_heading in mesh_headings:
            descriptor = mesh_heading.DescriptorName
            if descriptor and descriptor.get('MajorTopicYN') == 'Y':
                symptom_terms.add(descriptor.text)

    return symptom_terms

# Fetch symptom terms related to the disease
symptoms_related_to_disease = fetch_symptom_terms(disease)

# Print the symptom terms
print(f"Symptom terms related to the disease '{disease}':")
for symptom in symptoms_related_to_disease:
    print(symptom)

## to get top 100 abstracts from PubMed API
# from Bio import Entrez

# def get_symptoms_from_pubmed(disease_query):
#     """
#     Retrieves symptom-related information from PubMed based on a disease query.

#     Args:
#         disease_query (str): The disease name or keyword.

#     Returns:
#         list of str: List of symptom-related abstracts.
#     """
#     # Entrez.email = "your@email.com"  # Set your email (required by PubMed)

#     # Search for articles related to the disease
#     handle = Entrez.esearch(db="pubmed", term=disease_query, retmax=100)
#     record = Entrez.read(handle)
#     pmids = record["IdList"]

#     # Fetch abstracts for each article
#     symptom_abstracts = []
#     for pmid in pmids:
#         handle = Entrez.efetch(db="pubmed", id=pmid, retmode="text", rettype="abstract")
#         abstract = handle.read().strip()
#         symptom_abstracts.append(abstract)

#     return symptom_abstracts

# # Example usage
# disease_name = "COVID-19"
# symptom_abstracts = get_symptoms_from_pubmed(disease_name)

# # Print the symptom-related abstracts
# for i, abstract in enumerate(symptom_abstracts, start=1):
#     print(f"Abstract {i}:\n{abstract}\n")
