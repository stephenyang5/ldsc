import pandas as pd

df = pd.read_csv("uc_build37_45975_20161107.txt", sep="\t")
# Split MarkerName: '1:100000012_G_T' → chr=1, pos=100000012, ref=G, alt=T
variants = df["MarkerName"].str.extract(r"(?P<CHR>\d+):(?P<BP>\d+)_?(?P<REF>[ACGT])_?(?P<ALT>[ACGT])")
variants.to_csv("variants_for_lookup.txt", sep="\t", header=False, index=False)
