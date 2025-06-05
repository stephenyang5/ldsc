import pandas as pd

# Load your summary stats
df = pd.read_csv("cd_build37_40266_20161107.txt", sep="\t")

# Parse MarkerName into CHROM, POS, REF, ALT
marker_split = df['MarkerName'].str.extract(r'(?P<CHROM>[^:]+):(?P<POS>[0-9]+)_(?P<REF>[ACGT]+)_(?P<ALT>[ACGT]+)')

# Construct a valid VCF dataframe
df_vcf = pd.DataFrame({
    '#CHROM': marker_split['CHROM'],
    'POS': marker_split['POS'].astype(int),
    'ID': '.',
    'REF': marker_split['REF'],
    'ALT': marker_split['ALT'],
    'QUAL': '.',
    'FILTER': '.',
    'INFO': '.'
})

vcf_df = df_vcf[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

print(vcf_df.head())

# Write full VCF with header and contig lines
with open("cd_build37_40266_20161107_rsid.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for i in range(1, 23):
        f.write("##contig=<ID={0}>\n".format(i))
    f.write("##contig=<ID=X>\n")
    f.write("##contig=<ID=Y>\n")
    vcf_df.to_csv(f, sep="\t", index=False)

