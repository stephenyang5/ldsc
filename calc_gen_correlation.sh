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
--sumstats RA_GWASmeta_TransEthnic_v2.txt \
--N 58284 \
--out rart_okada_trans \
--snp SNPID \
--a1 A1 \
--a2 A2 \
--p P-val \
--signed-sumstats ORA1,1 \
--merge-alleles w_hm3.snplist


./munge_sumstats.py \
--out results/sumstats/diab1 \
--N 520580 \
--a1 effect_allele \
--a2 other_allele \
--snp variant_id \
--sumstats gwas_data/GCST90014023_buildGRCh38.tsv \
--merge-alleles w_hm3.snplist \
--p p_value

#arthritis and diabetes
python ldsc.py \
--rg rart_okada_european.sumstats.gz,output_file.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out rart_okada_diab1_eur_eur

#arthritis and diabetes
python ldsc.py \
--rg rart_okada_trans.sumstats.gz,output_file.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out rart_okada_diab1_trans_eur

python ldsc.py \
--rg rart_okada_trans.sumstats.gz,output_file.gz \
--ref-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--w-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--out rart_okada_diab1_trans

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



#Celiac Disease
./munge_sumstats.py \
--sumstats gwas_data/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv \
--N 15283 \
--out cel \
--snp rsid \
--a1 effect_allele \
--a2 other_allele \
--merge-alleles w_hm3.snplist \
--signed-sumstats OR,1 

python ldsc.py \
--rg sumstats/diab1.sumstats.gz,sumstats/cel.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out diab1_cel_eur

./munge_sumstats.py \
  --sumstats gwas_data/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv \
  --out results/sumstats/sle \
  --snp rsid \
  --a1 effect_allele \
  --a2 other_allele \
  --p p \
  --signed-sumstats OR,1 \
  --merge-alleles w_hm3.snplist \
  --N 18264

python ldsc.py \
--rg results/sumstats/uc.sumstats.gz,results/sumstats/uc.sumstats.gz \
--ref-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--w-ld-chr 1000GP_Phase3_LD/CEU_LD_scores.chr \
--out results/gen_correlation/uc_uc_eur

#ulcerative colitis
./munge_sumstats.py \
  --sumstats uc_gwas_with_rsid_test.tsv \
  --out results/uc_test \
  --snp rsID \
  --a1 allele2 \
  --a2 allele1 \
  --p P.value \
  --signed-sumstats Effect,0 \
  --N 18264

  #multiple sclerosis
  ./munge_sumstats.py \
  --sumstats gwas_data/imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv \
  --out results/sumstats/ms \
  --snp rsid \
  --a1 effect_allele \
  --a2 other_allele \
  --p p \
  --N 38582 \
  --ignore "OR"

#ulcerative colitis
  ./munge_sumstats.py \
  --sumstats gwas_data/uc_gwas_with_rsid.tsv \
  --out results/sumstats/uc \
  --snp rsid \
  --a1 Allele2 \
  --a2 Allele1 \
  --p P.value \
  --N 47281 \
  --ignore "OR,MarkerName"

  #pernicious anemia
  ./munge_sumstats.py \
  --sumstats gwas_data/pernicious_anemia_Laisketal2021_sumstats \
  --out results/sumstats/pa \
  --snp rs_number \
  --a1 reference_allele \
  --a2 other_allele \
  --frq eaf \
  --p p-value \
  --N-cas 2166\
  --N-con 659516 \
  --signed-sumstats z,0



#crohns disease sample size 47109
  ./munge_sumstats.py \
  --sumstats gwas_data/cd_gwas_with_rsid.tsv \
  --out results/sumstats/cd \
  --snp rsID \
  --a1 Allele2 \
  --a2 Allele1 \
  --p P.value \
  --N 47109 \
  --signed-sumstats Effect,0 \
  --ignore "OR,MarkerName"


sumstats_list=("diab1" "cd" "cel" "psoria" "rart" "sle" "uc" "ms" "pa")
for item in "${sumstats_list[@]}"; do
  python ldsc.py \
  --rg results/sumstats/${item}.sumstats.gz,results/sumstats/pa.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out results/gen_correlation/${item}_pa_eur
done