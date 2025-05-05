import csv

input_file = "GCST90472771.tsv"
pval_column = "p_value"

with open(input_file, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for i, row in enumerate(reader, start=2):  # start=2 to account for header
        pval = row.get(pval_column, "").strip()
        try:
            float(pval)
        except ValueError:
            print("Line {}: Invalid p-value: {}".format(i, pval))