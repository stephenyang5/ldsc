#List of current commands that I have been using to calculate genetic correlations between different diseases using LDSC.

#NEED TO WRITE THE WGET SCRIPTS FOR THIS AS WELL

# First downloading hapmap3 snplist to filter through SNPs for wrong imputation and standardization
wget https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist

# Clean sum stats files using hapmap3 snplist as reference
#for rheumatoid arthiritis
./munge_sumstats.py \
--sumstats GCST90132223_buildGRCh37.tsv \
--N 97173 \
--out rart \
--snp variant_id \
--a1 effect_allele \
--a2 other_allele \
--merge-alleles w_hm3.snplist

./munge_sumstats.py \
--out diab1 \
--N 501638 \
--a1 effect_allele \
--a2 other_allele \
--snp variant_id \
--sumstats GCST90014023_buildGRCh38.tsv \
--merge-alleles w_hm3.snplist \
--p p_value

#arthritis and diabetes
python ldsc.py \
--rg rart.sumstats.gz,output_file.gz \
--ref-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--w-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--out rart_diab1

./munge_sumstats.py \
--sumstats GCST90472771_rsid.tsv \
--N 494544 \
--out psoria \
--snp rsID \
--a1 EA \
--a2 NEA \
--merge-alleles w_hm3.snplist

#arthritis and psoriasis
python ldsc.py \
--rg rart.sumstats.gz,psoria.sumstats.gz \
--ref-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--w-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--out rart_psoria

#psoriasis and diabetes
python ldsc.py \
--rg output_file.gz,psoria.sumstats.gz \
--ref-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--w-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--out psoria_diab1

#download sumstats for crohns disease
curl -O ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2016-11-07/cd_build37_40266_20161107.txt.gz

#working with file to relabel
bgzip cd_build37_40266_20161107_rsid.vcf
tabix -p vcf cd_build37_40266_20161107_rsid.vcf.gz

#annotate the file with rsIDs usig markertxt_to_vcf.py

#zip file
bgzip cd_build37_40266_20161107_rsid.withcontigs.vcf

bcftools sort cd_build37_40266_20161107_rsid.vcf.gz \
    -Oz -o cd_build37_40266_20161107_rsid.sorted.vcf.gz

tabix -p vcf cd_build37_40266_20161107_rsid.sorted.vcf.gz

bcftools annotate \
  -a GCF_000001405.25.renamed.vcf.gz \
  -c CHROM,POS,ID,REF,ALT \
  cd_build37_40266_20161107_rsid.sorted.vcf.gz \
  -Oz -o cd_build37_40266_20161107_rsid.sorted.annotated.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' cd_build37_40266_20161107_rsid.sorted.vcf.gz > cd_build37_40266_20161107_extracted.tsv




./munge_sumstats.py \
--sumstats crohns_disease.tsv \
--N \
--out psoria \
--snp rsID \
--a1 EA \
--a2 NEA \
--merge-alleles w_hm3.snplist


