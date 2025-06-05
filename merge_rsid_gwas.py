import pandas as pd

# Load original GWAS
df = pd.read_csv("cd_build37_40266_20161107.txt", sep="\t")

# Parse MarkerName into CHROM, POS, REF, ALT
parts = df['MarkerName'].str.extract(r'(?P<CHROM>[^:]+):(?P<POS>[0-9]+)_(?P<REF>[ACGT]+)_(?P<ALT>[ACGT]+)', expand=True)
df['CHROM'] = parts['CHROM']
df['POS'] = parts['POS'].astype(int)
df['REF'] = parts['REF']
df['ALT'] = parts['ALT']

# Load rsID mapping from annotated VCF
rsids = pd.read_csv("rsid_map.tsv", sep="\t", header=None)
rsids.columns = ["CHROM", "POS", "REF", "ALT", "rsID"]

# Merge on 4 keys: CHROM, POS, REF, ALT
merged = pd.merge(df, rsids, how="left", on=["CHROM", "POS", "REF", "ALT"])

# Save to file
merged.to_csv("gwas_with_rsid.tsv", sep="\t", index=False)
