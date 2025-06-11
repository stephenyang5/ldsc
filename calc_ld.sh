#!/bin/bash
set -e

#create dir for the data
mkdir -p 1000G_CEU_data

#download vcf files for chromosomes 1–22
for chr in {1..22}; do
    wget -P 1000G_CEU_data ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    wget -P 1000G_CEU_data ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
done

#download dir
wget -P 1000G_CEU_data ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel


#filter for CEU data
awk '$2 == "CEU" {print $1}' 1000G_CEU_data/integrated_call_samples_v3.20130502.ALL.panel > 1000G_CEU_data/ceu_ids.txt

#Filter and process VCF files for CEU samples
mkdir -p 1000G_Phase3_VCFs

for chr in {1..22}; do
  bcftools view \
    -S 1000G_CEU_data/ceu_ids.txt \
    -O z \
    -o 1000G_Phase3_VCFs/CEU.chr${chr}.vcf.gz \
    1000G_Phase3_VCFs/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done

#convert to .bed/.bim/.fam format using plink2
for chr in {1..22}; do
        plink2 \
        --vcf eur_w/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --max-alleles 2 \
        --make-bed \
        --out eur_w/eur.chr${chr}
done


#download genetic map files for chromosomes 1–22
mkdir -p 1000G_Phase3

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/1000GP_Phase3.tgz
tar -xzf 1000GP_Phase3.tgz

wget -P 1000G_Phase3 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/1000GP_Phase3.tgz
tar -xzf 1000G_Phase3/1000GP_Phase3.tgz -C 1000G_Phase3

cd ..

mkdir -p 1000G_Phase3_CM

#annotate .bim files with cm map positions
for chr in {1..22}; do
    plink \
        --bfile 1000G_Phase3_VCFs/CEU.chr${chr} \
        --cm-map 1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt ${chr} \
        --make-bed \
        --out 1000G_Phase3_CM/CEU.chr${chr}.cm
done

#calculate LD scores using ldsc
#program is linear in snps and samples, and only uses 9 cores
# thus i pasted same command for different chromosomes into different terminal windows
for chr in {1..5}; do
    python ldsc.py \
        --bfile 1000G_Phase3_CM/CEU.chr${chr}.cm \
        --l2 \
        --ld-wind-cm 1 \
        --out CEU_LD_scores.chr${chr} \
        --yes-really \
        --thin-annot
done

for chr in {6..10}; do
    python ldsc.py \
        --bfile 1000G_Phase3_CM/CEU.chr${chr}.cm \
        --l2 \
        --ld-wind-cm 1 \
        --out CEU_LD_scores.chr${chr} \
        --yes-really \
        --thin-annot
done

for chr in {11..16}; do
    python ldsc.py \
        --bfile 1000G_Phase3_CM/CEU.chr${chr}.cm \
        --l2 \
        --ld-wind-cm 1 \
        --out CEU_LD_scores.chr${chr} \
        --yes-really \
        --thin-annot
done

for chr in {17..22}; do
    python ldsc.py \
        --bfile 1000G_Phase3_CM/CEU.chr${chr}.cm \
        --l2 \
        --ld-wind-cm 1 \
        --out CEU_LD_scores.chr${chr} \
        --yes-really \
        --thin-annot
done

#bcftools view 1000G_noCM_data/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '^#' | cut -f3 | head

for chr in {1..5}; do
    python ldsc.py \
        --bfile eur_w_CM/eur.chr${chr}.cm \
        --l2 \
        --ld-wind-cm 1 \
        --out eur_w_ld/chr${chr} \
        --yes-really \
        --thin-annot
done