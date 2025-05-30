#create dir for the data
mkdir -p 1000G_CEU_data
cd 1000G_CEU_data

#download vcf files for chromosomes 1–22
for chr in {1..22}; do
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
done

#download dir
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

#filter for CEU data
awk '$2 == "CEU" {print $1}' integrated_call_samples_v3.20130502.ALL.panel > ceu_ids.txt

#convert to .bed/.bim/.fam format using plink2
for chr in {1..22}; do
        plink2 \
        --vcf ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --keep ceu_ids.txt \
        --max-alleles 2\
        --make-bed \
        --out CEU.chr${chr}
done

#fill in cM map data
for chr in {1..22}; do
    plink \
        --bfile 1000G_CEU_data/CEU.chr${chr} \
        --cm-map 1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt ${chr} \
        --make-bed \
        --out CEU.chr${chr}.cm
done


#calculate LD scores using ldsc
#program is linear in snps and samples, and only uses 9 cores
# thus i pasted same command for different chromosomes into different terminal windows
for chr in {13..22}; do
    python ldsc.py \
        --bfile 1000G_CEU_CM_data/CEU.chr${chr}.cm \
        --l2 \
        --ld-wind-cm 1 \
        --out CEU_LD_scores.chr${chr} \
        --yes-really \
        --thin-annot
done

