import os
import json
import faiss
import numpy as np

from icd9_vector_db import model, icd9_index, icd9_terms

# Function to load the FAISS index
def load_faiss_index(index_path):
    return faiss.read_index(index_path)

def search_icd9(query, k=1):
    query_embedding = model.encode([query]).astype('float32')
    D, I = icd9_index.search(query_embedding, k)
    neighbors_list=[]
    
    for i in range(len(I[0])):
        res=icd9_terms[I[0][i]]
        res+=" - " + str(D[0][i])
        neighbors_list.append(res)
    return neighbors_list

# Function to get the nearest term from the FAISS index
def get_nearest_term(faiss_index, term, term_to_vector_mapping):
    if term not in term_to_vector_mapping:
        raise ValueError(f"Term '{term}' not found in the term-to-vector mapping.")
    vector = term_to_vector_mapping[term]
    vector = np.array([vector], dtype=np.float32)
    D, I = faiss_index.search(vector, 1)  # Search for the nearest neighbor
    nearest_term_index = I[0][0]
    return nearest_term_index

# Directory containing the relation directories
relation_dir = 'relationships'

# Initialize a dictionary to hold the results
results = {}

# Traverse the directories and process relation_count.txt files
for disease_dir in os.listdir(relation_dir):
    disease_dir_path = os.path.join(relation_dir, disease_dir)
    print(disease_dir_path)
    if os.path.isdir(disease_dir_path):
        relation_count_file = os.path.join(disease_dir_path, f'{disease_dir}_relation_count.txt')
        term = " ".join(disease_dir.split('_'))
        if os.path.exists(relation_count_file):
            with open(relation_count_file, 'r') as f:
                for line in f:
                    mapping_list = line.strip().split(', ')
                    neg_weight = mapping_list[-1]
                    pos_weight = mapping_list[-2]
                    if len(mapping_list)>4:
                        if len(disease_dir.split('_'))>1:
                            node1=", ".join(mapping_list[0:len(disease_dir.split('_'))])
                            node2 = ", ".join(mapping_list[len(disease_dir.split('_')):-2])
                        else:
                            node1 = mapping_list[0]
                            node2 = ", ".join(mapping_list[1:-2])
                    else:
                        node1 = mapping_list[-4]
                        node2 = mapping_list[-3]
                    nearest_term = search_icd9(node2,k=1)
                    if node1 not in results:
                        results[node1] = []
                    # Check if the nearest_term is already mapped to node1
                    existing_term = next((item for item in results[node1] if item['disease'] == nearest_term), None)
                    
                    if existing_term:
                        if existing_term['pos_weight'] > pos_weight:
                            existing_term['pos_weight'] = pos_weight
                    else:
                        results[node1].append({'disease': ''.join(nearest_term), 'pos_weight': pos_weight})

# Save the results to a JSON file
output_file = 'mapped_terms_icd9.json'
with open(output_file, 'a') as f:
    json.dump(results, f, indent=4)

print(f'Results saved to {output_file}')