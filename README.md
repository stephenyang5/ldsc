# Personal Work
This is documentation for exploratory work looking at genetic correlation between AIDs. 

All the data is publicly available to download. If you cannot find the 'wget' function call or link that I used to download a piece of data let me know and I will update it!

# Formatting Summary Statistic Data
This repo is forked from LDSC, a program written for various operation related to LD scores. It relies on inputted GWAS summary statistics having an rsID identifier in addition to of chromosome and position locations. However, not all raw summary statistic files contained this column and I therefore needed to go through and clean some datasets. The following is my workflow. 

## Crohn's Disease and Ulcerative Colitis
The summary statistics for these two files contain Marker identifiers instead of rsID, chromosome, and position. These are given in the format: 1:100000012_G_T
which corresponds to chromosome, position, effect allele, other allele. 

Therefore, I ran a regex over this column to parse out the required portions, queried them against an rsID reference file, and then merged the queried results into the raw data (producing a new output file). The methods for this can be found in 'marker_cols_conversion.sh'.

Note that the rsID reference file gives chromosome numbers in RefSeq format e.g. NC_000001.10	is the name of chromosome 1. The process for this is also given in the 'marker_cols_conversion.sh' file. (THIS NEEDS TO BE DONE STILL)


# Calculating LD scores 
As an Excercise I Calculated LD scores in the CEU sub population of the 1000G EUR population. All of the code to do this can be seen in the calc_ld.sh with the exact commands that I used to operate on the data. While I was working, I had a /data/ subfolder in my working directory, but I removed folder organization from these scripts and it is not standardized. I will go back and fix so that once all data is downloaded, then the script can be run start to finish.

For the summary statistics I downloaded, please see the following link to a spreadsheet that has all the data that I used

# Calculating Genetic Correlation and LD Scores

For the graphs generated I used pre-computed LD scores on the entire EUR population from the creators of LDSC. Their link to download is since deprecated. I downloaded these LD scores from [this location](https://zenodo.org/records/8182036).

## 'calc_gen_correlation.sh'
This file contains all of the bash calls that I used in preparing my summary statistics for LDSC analysis. 

'munge_sumstats.py' is the build in ldsc function used to standardize summary statistics to their formatting. There could be a clever way to loop through and munge, but since column names such as effect_allele, and population size vary between studies, I wrote this out manually. There is a loop for the actual genetic correlation calculation at the very end once data has been cleaned.

## plot_heatmap.py
After working through the genetic correlation file, this plotting script will visualize the results in a heatmap.

 
# For Sumstats in the form 1:100000012_G_T	

First use an AWK command to run a regex and convert into CHROM ID SNP format.

e.g.

# Extract CHR, POS from your file
awk 'BEGIN{FS=OFS="\t"} NR>1 {
  split($1,a,":|_");
  print a[1], a[2], a[3], a[4];
}' gwas_data/.txt > cd_build37_40266_20161107.txt

Side note, you need to annotate GCF_000001405 with chromosome IDs since it is in RefSeq format



Then query this file to collate rsIDs  GCST has chromosome names in the NC refseq format. Therefore, we must rename them to numerical chromosome numbers which is what bcft tools recognizes. After doing this, can query the GCST reference file for rsids.

This command renames the chromosomes
bcftools annotate --rename-chrs chr_rename.txt -O z -o renamed.GCF_000001405.25.vcf.gz GCF_000001405.25.gz

Query renamed gcf file for rsids of the chrom/pos found in GWAS
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' -R markers_for_query.tsv All_20180418.vcf.gz > rsid_lookup.tsv



Finally, use pandas operations to insert rsIDs into original file using merge_rsid_gwas.txt

# LDSC (LD SCore) `v1.0.1`

`ldsc` is a command line tool for estimating heritability and genetic correlation from GWAS summary statistics. `ldsc` also computes LD Scores.

## Getting Started



In order to download `ldsc`, you should clone this repository via the commands
```  
git clone https://github.com/bulik/ldsc.git
cd ldsc
```

In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with LDSC's dependencies:

```
conda env create --file environment.yml
source activate ldsc
```

Once the above has completed, you can run:

```
./ldsc.py -h
./munge_sumstats.py -h
```
to print a list of all command-line options. If these commands fail with an error, then something as gone wrong during the installation process. 

Short tutorials describing the four basic functions of `ldsc` (estimating LD Scores, h2 and partitioned h2, genetic correlation, the LD Score regression intercept) can be found in the wiki. If you would like to run the tests, please see the wiki.

## Updating LDSC

You can update to the newest version of `ldsc` using `git`. First, navigate to your `ldsc/` directory (e.g., `cd ldsc`), then run
```
git pull
```
If `ldsc` is up to date, you will see 
```
Already up-to-date.
```
otherwise, you will see `git` output similar to 
```
remote: Counting objects: 3, done.
remote: Compressing objects: 100% (3/3), done.
remote: Total 3 (delta 0), reused 0 (delta 0), pack-reused 0
Unpacking objects: 100% (3/3), done.
From https://github.com/bulik/ldsc
   95f4db3..a6a6b18  master     -> origin/master
Updating 95f4db3..a6a6b18
Fast-forward
 README.md | 15 +++++++++++++++
 1 file changed, 15 insertions(+)
 ```
which tells you which files were changed. If you have modified the `ldsc` source code, `git pull` may fail with an error such as `error: Your local changes to the following files would be overwritten by merge:`. 

In case the Python dependencies have changed, you can update the LDSC environment with

```
conda env update --file environment.yml
```

## Where Can I Get LD Scores?

You can download [European](https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2) and [East Asian LD Scores](https://data.broadinstitute.org/alkesgroup/LDSCORE/eas_ldscores.tar.bz2) from 1000 Genomes [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/). These LD Scores are suitable for basic LD Score analyses (the LD Score regression intercept, heritability, genetic correlation, cross-sex genetic correlation). You can download partitioned LD Scores for partitioned heritability estimation [here](http://data.broadinstitute.org/alkesgroup/LDSCORE/).


## Support

Before contacting us, please try the following:

1. The [wiki](https://github.com/bulik/ldsc/wiki) has tutorials on [estimating LD Score](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial), [heritability, genetic correlation and the LD Score regression intercept](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) and [partitioned heritability](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability).
2. Common issues are described in the [FAQ](https://github.com/bulik/ldsc/wiki/FAQ)
2. The methods are described in the papers (citations below)

If that doesn't work, you can get in touch with us via the [google group](https://groups.google.com/forum/?hl=en#!forum/ldsc_users).

Issues with LD Hub?  Email ld-hub@bristol.ac.uk


## Citation

If you use the software or the LD Score regression intercept, please cite

[Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)

For genetic correlation, please also cite

[Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3406) Preprint available on bioRxiv doi: http://dx.doi.org/10.1101/014498

For partitioned heritability, please also cite

[Finucane, HK, et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3404) Preprint available on bioRxiv doi: http://dx.doi.org/10.1101/014241

For stratified heritability using continuous annotation, please also cite

[Gazal, S, et al. Linkage disequilibrium–dependent architecture of human complex traits shows action of negative selection. Nature Genetics, 2017.](https://www.nature.com/articles/ng.3954) 

If you find the fact that LD Score regression approximates HE regression to be conceptually useful, please cite

Bulik-Sullivan, Brendan. Relationship between LD Score and Haseman-Elston, bioRxiv doi: http://dx.doi.org/10.1101/018283

For LD Hub, please cite

[Zheng, et al. LD Hub: a centralized database and web interface to perform LD score regression that maximizes the potential of summary level GWAS data for SNP heritability and genetic correlation analysis. Bioinformatics (2016)](https://doi.org/10.1093/bioinformatics/btw613)


## License

This project is licensed under GNU GPL v3.


## Authors

Brendan Bulik-Sullivan (Broad Institute of MIT and Harvard)

Hilary Finucane (MIT Department of Mathematics)
