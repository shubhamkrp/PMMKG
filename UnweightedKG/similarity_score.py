def levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def similarity_score(s1, s2):
    lev_distance = levenshtein_distance(s1, s2)
    max_len = max(len(s1), len(s2))
    score = (1 - lev_distance / max_len) * 100
    return score

term1 = "depression"
term2 = "chord compression"
score = similarity_score(term1, term2)

print(f"Similarity score between '{term1}' and '{term2}': {score:.2f}%")
