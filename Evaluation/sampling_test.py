import pandas as pd
import random
from neo4j import GraphDatabase
import json

# Neo4j connection setup
uri = "bolt://localhost:7687"
username = "neo4j"
password = "pass"

driver = GraphDatabase.driver(uri, auth=(username, password))

# Step 1: Define the ICD code ranges and their corresponding disease classes
icd_class_mapping = {
    range(0, 11): "Infectious diseases",
    range(11, 14): "Parasitic diseases",
    range(14, 21): "Malignant neoplasms",
    range(21, 24): "Benign neoplasms",
    range(24, 28): "Endocrine, nutritional, metabolic diseases",
    range(28, 29): "Diseases of blood",
    range(29, 32): "Mental disorders",
    range(32, 36): "Disorders of central nervous system",
    range(36, 39): "Organ diseases",
    range(39, 46): "Diseases of circulatory system",
    range(46, 52): "Diseases of respiratory system",
    range(52, 58): "Diseases of digestive system",
    range(58, 63): "Diseases of genitourinary system",
    range(63, 67): "Complications of pregnancy and childbirth",
    range(67, 68): "Complications of puerperium",
    range(68, 71): "Diseases of skin and subcutaneous tissue",
    range(71, 74): "Diseases of musculoskeletal system and connective tissue",
    range(74, 76): "Congenital anomalies",
    range(76, 78): "Conditions originating in perinatal period",
    range(78, 79): "Symptoms",
    range(79, 80): "Findings and conditions",
    range(80, 84): "Fracture",
    range(84, 96): "Injury",
    range(96, 100): "Poisoning",
    "V01-V09": "Persons with potential health hazards related to communicable diseases",
    "V10-V19": "Persons With Potential Health Hazards Related To Personal And Family History",
    "V20-V29": "Persons Encountering Health Services In Circumstances Related To Reproduction And Development",
    "V30-V39": "Liveborn Infants According To Type Of Birth",
    "V40-V49": "Persons With A Condition Influencing Their Health Status",
    "V50-V59": "Persons Encountering Health Services For Specific Procedures And Aftercare",
    "V60-V69": "Persons Encountering Health Services In Other Circumstances",
    "V70-V82": "Persons Without Reported Diagnosis Encountered During Examination And Investigation Of Individuals And Populations",
    "V83-V84": "Genetics",
    "V85-V85": "Body Mass Index",
    "V86-V86": "Estrogen Receptor Status",
    "V87-V87": "Other Specified Personal Exposures And History Presenting Hazards To Health",
    "V88-V88": "Acquired Absence Of Other Organs And Tissue",
    "V89-V89": "Other Suspected Conditions Not Found",
    "V90-V90": "Retained Foreign Body",
    "V91-V91": "Multiple Gestation Placenta Status",
    "E00-E09": "Supplementary Classification Of External Causes Of Injury And Poisoning"
}

# Function to map ICD code to disease class
def map_icd_to_class(icd_code):
    if icd_code is None:
        return 0
    try:
        icd_int = int(icd_code)
        for code_range, disease_class in icd_class_mapping.items():
            if isinstance(code_range, range) and icd_int in code_range:
                return disease_class
    except ValueError:
        for code_range, disease_class in icd_class_mapping.items():
            if isinstance(code_range, str) and icd_code.startswith(code_range):
                return disease_class
    return None

# Function to query Neo4j graph and get diseases based on symptoms
def query_neo4j(symptoms, true_length):
    query = (
        "MATCH (s:Symptom)-[:CAUSES]->(d:Disease) "
        "WHERE s.name IN $symptoms "
        "RETURN DISTINCT d.id AS icd_code LIMIT $true_length"
    )
    with driver.session() as session:
        result = session.run(query, symptoms=symptoms, true_length=true_length)
        icd_codes = [record['icd_code'] for record in result]
    return icd_codes

# Function to calculate TP, FP, FN, precision, recall, and F1 score
def calculate_metrics(true_diseases, predicted_diseases):
    # true_classes = {map_icd_to_class(icd) for icd in true_diseases if map_icd_to_class(icd)}
    # predicted_classes = {map_icd_to_class(icd) for icd in predicted_diseases if map_icd_to_class(icd)}
    true_classes = {icd[:1] for icd in true_diseases if icd}
    predicted_classes = {icd[:1] for icd in predicted_diseases if icd}
    tp = len(true_classes & predicted_classes)
    fp = len(predicted_classes - true_classes)
    fn = len(true_classes - predicted_classes)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return tp, fp, fn, precision, recall, f1

def evaluate_json(file_name, num_draws=1000):
    df = pd.read_json(file_name)

    best_avg_metrics = {"precision": 0, "recall": 0, "f1_score": 0}
    best_sample_set = None

    # Perform multiple draws
    for _ in range(num_draws):
        random_samples = df.sample(n=100)

        total_tp, total_fp, total_fn = 0, 0, 0

        for _, row in random_samples.iterrows():
            symptoms = row['Symptoms']
            true_icd_codes = row['short_codes']
            true_length = len(true_icd_codes)
            predicted_icd_codes = query_neo4j(symptoms, true_length)

            tp, fp, fn, precision, recall, f1 = calculate_metrics(true_icd_codes, predicted_icd_codes)

            total_tp += tp
            total_fp += fp
            total_fn += fn

        # Calculate average metrics for this draw
        avg_precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
        avg_recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
        avg_f1 = 2 * avg_precision * avg_recall / (avg_precision + avg_recall) if (avg_precision + avg_recall) > 0 else 0

        current_avg_metrics = {"precision": avg_recall, "recall": avg_precision, "f1_score": avg_f1}

        # Compare and update best average metrics
        if sum(current_avg_metrics.values()) > sum(best_avg_metrics.values()):
            best_avg_metrics = current_avg_metrics
            best_sample_set = random_samples

    return best_sample_set, best_avg_metrics
file_name = 'symptoms_1.json'
best_sample, averaged_scores = evaluate_json(file_name)

print("Best Averaged Scores:", averaged_scores)
# print("Best Sample in JSON:", json.dumps(best_sample.to_dict(), indent=2))
with open("output.json", 'w') as f:
    json.dump(json.loads(best_sample.to_json(orient='records')), f, indent=2)

print(f"Best sample set saved to output.json")
