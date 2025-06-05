import pandas as pd
import gwaslab as gl

def main():
    # Set up gwaslab references
    gl.download_ref("1kg_dbsnp151_hg19_auto")
    ref_tsv = gl.get_path("1kg_dbsnp151_hg19_auto")
    ref_vcf = "GCF_000001405.25.renamed.vcf.gz"
    chr_dict = gl.get_number_to_NC(build="19")

    # Stream and preprocess the summary stats
    chunk_iter = pd.read_csv("cd_build37_40266_20161107.txt", sep="\t", chunksize=100000)

    for i, chunk in enumerate(chunk_iter):
        print("Processing chunk {}".format(i))

        # Parse MarkerName into CHR, POS, REF, ALT
        parts = chunk["MarkerName"].str.extract(r'([^:]+):([0-9]+)_([ACGT]+)_([ACGT]+)', expand=True)
        chunk["chromosome"] = parts[0]
        chunk["base_pair_location"] = parts[1].astype(int)
        chunk["other_allele"] = parts[2]  # REF = other
        chunk["effect_allele"] = parts[3]  # ALT = effect

        # Create GWAS object
        mysumstats = gl.Sumstats(
            chunk,
            chrom="chromosome",
            pos="base_pair_location",
            ea="effect_allele",
            nea="other_allele",
            beta="Effect",
            se="StdErr",
            p="P.value",
            direction="Direction"
        )

        # Clean, harmonize, and map rsIDs
        mysumstats.basic_check()
        mysumstats.fix_id(fixsep=True)
        mysumstats.assign_rsid(
            ref_rsid_tsv=ref_tsv,
            ref_rsid_vcf=ref_vcf,
            chr_dict=chr_dict,
            n_cores=4
        )

        # Output chunk to combined file
        mysumstats.data.to_csv("cd_build37_40266_20161107_rsid.tsv", sep="\t", index=False,
                               mode='a' if i > 0 else 'w', header=(i == 0))

if __name__ == "__main__":
    main()