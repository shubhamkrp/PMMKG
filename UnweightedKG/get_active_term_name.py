import csv
import argparse
import os

def extract_name_from_terms(fname):
    data = open(fname, 'r', encoding="utf8")
    # print(data)

    # Split data into individual terms
    terms = data.read().split("\n\n")
    # print(len(terms))

    # Initialize a list to store extracted names
    extracted_names = []

    # Process each term
    for term in terms:
        lines = term.split("\n")
        name = None
        for line in lines:
            if line.startswith("name:"):
                name = line.split(":")[1].strip()
            elif line.startswith("is_a:"):
                extracted_names.append(name)
                break

    return extracted_names

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Parse .obo file to find term name.')
    parser.add_argument('input_file', help='Input file path. File should be the *.obo file.')
    parser.add_argument('--output_dir', help='Directory to output files', default='.')
    parser.add_argument('output_file', help='Output file path. File should be in csv format.')
    args = parser.parse_args()

    output_fname = os.path.join(args.output_dir, args.output_file)

    extracted_names = extract_name_from_terms(args.input_file)
    if len(extracted_names)==0:
        raise Exception("An exception occurred since parsed entity is empty")
    
    # Save extracted names to a CSV file
    with open(output_fname, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Name"])
        for name in extracted_names:
            writer.writerow([name])

    print(f"Extracted term names saved to {output_fname}")


#########---------USAGE-----------#############
#get_active_term_name.py [-h] [--output_dir OUTPUT_DIR] input_file output_file
###############################################