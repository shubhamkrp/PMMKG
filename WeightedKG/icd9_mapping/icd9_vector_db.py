import numpy as np
from sentence_transformers import SentenceTransformer
import faiss
import re

def read_icd9_cm_terms(file_path):
    icd9_terms = []
    with open(file_path, 'r', encoding="ISO-8859-1") as file:
        for line in file:
            parts = re.split(r'\s{2,}', line.strip())
            if len(parts) == 2:
                code, name = parts
                icd9_terms.append(f"{code} - {name}")
    return icd9_terms

file_path = 'icd9_cm_codes.txt'
icd9_terms = read_icd9_cm_terms(file_path)


model = SentenceTransformer('dmis-lab/biobert-base-cased-v1.1')

# Generate embeddings for ICD-9-CM terms
icd9_embeddings = model.encode(icd9_terms)

# Set Up FAISS Vector Database
# Convert embeddings to a format suitable for FAISS
icd9_embeddings = np.array(icd9_embeddings).astype('float32')

# Create a FAISS index
icd9_index = faiss.IndexFlatL2(icd9_embeddings.shape[1])
icd9_index.add(icd9_embeddings)

faiss.write_index(icd9_index, 'faiss_index_icd9.bin')
