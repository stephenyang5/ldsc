import gwaslab as gl
import pandas as pd

def main():
    gl.download_ref("1kg_dbsnp151_hg19_auto")
    ref_tsv = gl.get_path("1kg_dbsnp151_hg19_auto")
    ref_vcf = "GCF_000001405.25.renamed.vcf.gz"
    chr_dict = gl.get_number_to_NC(build="19")

    chunk_iter = pd.read_csv("GCST90472771.tsv", sep="\t", chunksize=100000)

    for i, chunk in enumerate(chunk_iter):
        print(f"Processing chunk {i}")
        mysumstats = gl.Sumstats(chunk, chrom="chromosome", pos="base_pair_location",
                                 ea="effect_allele", nea="other_allele", beta="beta",
                                 se="standard_error", p="p_value", direction="direction",
                                 n="cum_eff_sample_size")
        mysumstats.basic_check()
        mysumstats.fix_id(fixsep=True)
        mysumstats.assign_rsid(
            ref_rsid_tsv=ref_tsv,
            ref_rsid_vcf=ref_vcf,
            chr_dict=chr_dict,
            n_cores=8  # Avoid oversaturating memory/CPU
        )
        mysumstats.data.to_csv("GCST90472771_rsid.tsv", sep="\t", index=False,
                               mode='a' if i > 0 else 'w', header=(i == 0))

if __name__ == "__main__":
    main()