import csv
import re
import os
from collections import defaultdict
import spacy
from negspacy.negation import Negex
from spacy.pipeline import Sentencizer
from new_mining import term_name

# Function to read terms from the input text file
def read_terms(filename):
    with open(filename, 'r') as file:
        terms = [line.strip() for line in file.readlines()]
    return terms

# Function to create pairs of the first term with all the other terms
def create_pairs(terms):
    pairs = [(terms[0], terms[i]) for i in range(1, len(terms))]
    return pairs

# Function to read titles and abstracts from the CSV file
def read_articles(filename):
    articles = []
    with open(filename, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            articles.append({
                'title': row['title'],
                'abstract': row['abstract']
            })
    return articles

# Function to detect if a pair is negatively related in a sentence using negex
def is_negative(sentence, term1, term2):
    doc = nlp(sentence)
    for entity in doc.ents:
        if entity.text.lower() == term1.lower() or entity.text.lower() == term2.lower():
            if entity._.negex:
                return True
    return False

# Function to search for pairs in titles and abstracts sentence-wise
def count_co_occurrences(pairs, articles):
    pair_counts = defaultdict(lambda: {'positive_count': 0, 'negative_count': 0})

    for article in articles:
        text = f"{article['title']} {article['abstract']}"
        sentences = re.split(r'(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?)\s', text)

        for sentence in sentences:
            for pair in pairs:
                term1, term2 = pair
                if re.search(r'\b' + re.escape(term1) + r'\b', sentence) and re.search(r'\b' + re.escape(term2) + r'\b', sentence):
                    if is_negative(sentence, term1, term2):
                        pair_counts[pair]['negative_count'] += 1
                    else:
                        pair_counts[pair]['positive_count'] += 1

    return pair_counts

# Initialize spaCy model with negex
nlp = spacy.load("en_core_web_sm")
nlp.add_pipe("negex", config={"ent_types": ["NOUN", "PROPN", "VERB", "ADJ"]})
nlp.add_pipe("sentencizer")

if not os.path.exists(term_name):
    os.mkdir(term_name)
terms_filename = os.path.join(term_name, f'{term_name}.txt')
csv_filename = os.path.join(term_name, f'{term_name}_articles.csv')

# Read terms and create pairs
terms = read_terms(terms_filename)
pairs = create_pairs(terms)

# Read articles from CSV
articles = read_articles(csv_filename)

# Count co-occurrences of pairs
pair_counts = count_co_occurrences(pairs, articles)

with open(os.path.join(term_name, f'{term_name}_relation_count.txt'), 'w') as count_file:
    for pair, counts in pair_counts.items():
        count_file.write(f"{pair[0]}, {pair[1]}, {counts['positive_count']}, {counts['negative_count']}\n")
