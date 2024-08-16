import requests
from bs4 import BeautifulSoup
import pandas as pd
from similarity_bw_pubmed_do import similar


# Define the base URL for the PubMed API
PUBMED_API_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Define the function to fetch disease terms related to the given symptom
def fetch_disease_terms(symptom):
    # Replace spaces with '+' for URL encoding
    symptom_query = symptom.replace(" ", "+")

    # Construct the search query URL
    search_url = f"{PUBMED_API_BASE_URL}esearch.fcgi?db=pubmed&term={symptom_query}[MeSH Terms]&retmax=50"

    # Send the request to PubMed API
    response = requests.get(search_url)
    response.raise_for_status()  # Raise an error if the request failed

    soup = BeautifulSoup(response.content, 'xml')

    # Extract the PubMed IDs (PMIDs) from the response
    pmids = [id_tag.text for id_tag in soup.find_all('Id')]

    # Initialize a set to store unique disease terms
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


symptom_file="OUTPUT/active_symp_terms.csv"
disease_file="OUTPUT/active_disease_terms.csv"

symptom_df=pd.read_csv(symptom_file)
disease_df=pd.read_csv(disease_file, encoding="cp1252")

symptom_df=symptom_df.loc[80:100] #only to cross-check, remove this later on

#select name column
column_name="Name"
symptom_column=symptom_df[column_name]
disease_column=disease_df[column_name]

COLUMN_NAMES=["Symptom", "Disease"]
df=pd.DataFrame(columns=COLUMN_NAMES)

#traverse through each term
for symptom in symptom_column:

    # Fetch disease terms related to the symptom
    diseases_related_to_symptom = fetch_disease_terms(symptom)
    print(len(diseases_related_to_symptom))
    total_diseases = set(disease_column)

    actual_diseases=set()
    for d in diseases_related_to_symptom:
        for d_vocab in total_diseases:
            if similar(d,d_vocab):
                actual_diseases.add(d_vocab)
    # actual_diseases = list(diseases_related_to_symptom & total_diseases)
    actual_diseases_string = ';'.join(actual_diseases)
    df.loc[len(df)]=[symptom,actual_diseases_string]
    # print(len(actual_diseases))
print(df.head())


df.to_csv("OUTPUT/relationship3.csv",sep=",",mode='a', index=False, header=False)


    # # Print the disease terms
    # print(f"Disease terms related to the symptom '{symptom}':")
    # for disease in diseases_related_to_symptom:
    #     print(disease)
