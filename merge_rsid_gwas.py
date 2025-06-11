import pandas as pd

# Load original GWAS
df = pd.read_csv("gwas_data/cd_build37_40266_20161107.txt", sep="\t")

# Parse MarkerName into CHROM, POS, REF, ALT
parts = df['MarkerName'].str.extract(r'(?P<CHROM>[^:]+):(?P<POS>[0-9]+)_(?P<REF>[ACGT]+)_(?P<ALT>[ACGT]+)', expand=True)


df['CHROM'] = parts['CHROM'].astype(str)
df['POS'] = parts['POS'].astype(str)
df['REF'] = parts['REF'].astype(str)
df['ALT'] = parts['ALT'].astype(str)

# Load rsID mapping from annotated VCF
rsids = pd.read_csv("gwas_data/cd_rsid_lookup.tsv", sep="\t", header=None)
rsids.columns = ["CHROM", "POS", "REF", "ALT", "rsID"]

rsids['CHROM'] = rsids['CHROM'].astype(str)
rsids['POS'] = rsids['POS'].astype(str)
rsids['REF'] = rsids['REF'].astype(str)
rsids['ALT'] = rsids['ALT'].astype(str)

df['REF'] = df['REF'].str.upper()
df['ALT'] = df['ALT'].str.upper()
rsids['REF'] = rsids['REF'].str.upper()
rsids['ALT'] = rsids['ALT'].str.upper()

# Ensure CHROM format matches (e.g. '1' vs 'chr1')
rsids['CHROM'] = rsids['CHROM'].astype(str).str.replace('^chr', '', regex=True)
df['CHROM'] = df['CHROM'].astype(str)

# Make sure POS is string
df['POS'] = df['POS'].astype(str)
rsids['POS'] = rsids['POS'].astype(str)


# Merge on 4 keys: CHROM, POS, REF, ALT
merged = pd.merge(df, rsids, how="left", on=["CHROM", "POS", "REF", "ALT"])

# Save to file
merged.to_csv("gwas_data/cd_gwas_with_rsid.tsv", sep="\t", index=False)
