import json
import pandas as pd
from neo4j import GraphDatabase
from sklearn.metrics import precision_score, recall_score, f1_score

# Step 1: Load JSON file
with open("symptom_disease.json", "r") as file:
    data = json.load(file)

# Step 2: Setup Neo4j driver and session
uri = "bolt://localhost:7687"  # Modify according to your setup
driver = GraphDatabase.driver(uri, auth=("neo4j", "password"))

# Function to query Neo4j for diseases based on symptoms
def get_predicted_diseases(symptoms):
    query = """
    MATCH (s:Symptom)-[:CAUSES]->(d:Disease)
    WHERE s.name IN $symptoms
    RETURN DISTINCT d.name AS disease
    """
    with driver.session() as session:
        result = session.run(query, symptoms=symptoms)
        return [record["disease"][:5] for record in result]

# Function to check if the predicted disease falls within the range
def is_correct_classification(predicted, actual_ranges):
    pred_prefix = predicted[:2]
    for range_str in actual_ranges:
        if "-" in range_str:
            start, end = range_str.split("-")
            if start[:2] <= pred_prefix <= end[:2]:
                return True
        else:  # Single disease code case
            if pred_prefix == range_str[:2]:
                return True
    return False

# Step 3: Iterate through JSON data and calculate metrics
precisions = []
recalls = []
f1_scores = []

for entry in data:
    symptoms = entry["symptoms"]
    actual_ranges = entry["diseases"]
    
    predicted_diseases = get_predicted_diseases(symptoms)
    
    y_true = []
    y_pred = []

    for predicted in predicted_diseases:
        if is_correct_classification(predicted, actual_ranges):
            y_true.append(1)
            y_pred.append(1)
        else:
            y_true.append(0)
            y_pred.append(1)

    for actual_range in actual_ranges:
        if not any(is_correct_classification(pred, [actual_range]) for pred in predicted_diseases):
            y_true.append(1)
            y_pred.append(0)
    
    if y_true:
        precisions.append(precision_score(y_true, y_pred, zero_division=0))
        recalls.append(recall_score(y_true, y_pred, zero_division=0))
        f1_scores.append(f1_score(y_true, y_pred, zero_division=0))

# Step 4: Calculate Average Precision, Recall, and F1-score
avg_precision = sum(precisions) / len(precisions) if precisions else 0
avg_recall = sum(recalls) / len(recalls) if recalls else 0
avg_f1 = sum(f1_scores) / len(f1_scores) if f1_scores else 0

print(f"Average Precision: {avg_precision}")
print(f"Average Recall: {avg_recall}")
print(f"Average F1 Score: {avg_f1}")

# Step 5: Close Neo4j connection
driver.close()
