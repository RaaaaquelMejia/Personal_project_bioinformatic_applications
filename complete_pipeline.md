#######################################################
#    Inbreeding rate estimation of Numenius arquata   #
#######################################################

# Author : Raquel Mejia
# Last modification : 27/05/2026

#The main goal of this projects is to quantify the inbreeding
 of one endangered bird species called Numenius arquata.

## Quality control of data with genomepanel_nf ##

# Allows to install the pipeline
 module load genomepanel_nf

# Downloads directly the data of the reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/106/895/GCF_964106895.1_bNumArq3.hap1.1/GCF_964106895.1_bNumArq3.hap1.1_genomic.fna.gz
gunzip GCF_964106895.1_bNumArq3.hap1.1_genomic.fna.gz

# Subsetting reference genome to chromosomal contigs (starting with ">NC_...")
module load SAMtools

samtools faidx GCF_964106895.1_bNumArq3.hap1.1_genomic.fna
grep "^>NC_" GCF_964106895.1_bNumArq3.hap1.1_genomic.fna | awk '{print substr($1,2)}' | xargs samtools faidx GCF_964106895.1_bNumArq3.hap1.1_genomic.fna > GCF_964106895.1_bNumArq3.hap1.1_genomic.subset.fna

# Run of the pipeline
# Need of SRA.txt that includes the SRA of 54 Numenius arquata individuals
screen -S gp

  nextflow run $GENOMEPANEL_HOME/pipeline/main.nf \
     --reference /home/ba-student2/Personal_project_bioinformatic_applications/GCF_964106895.1_bNumArq3.hap1.1_genomic.subset.fna \
     --reference_segments 20000000 \
     --SRA_index /home/ba-student2/Personal_project_bioinformatic_applications/SRA.txt \
     --ploidy 2 \
     --outdir results_numenius -profile slurm

screen -r gp


## Filtering obtained data with VCFtools ##

# Downloading VCFtools 

module spider vcftools # to find hiding program already installed
module load VCFtools # load the programm

# VCFtools filtering with --keep option to use selected samples
# First, filter with --max-missing 0.8 (keep variants with 20% or less missing data)
# Using --stdout to pipe directly to bgzip for compression
vcftools --gzvcf results_numenius/final_variants.clean.vcf.gz --keep results_numenius/Selected_SRR.tsv --max-missing 0.8 --recode --recode-INFO-all --stdout | bgzip -c > results_numenius/final_variants.clean.max_missing_0.8.recode.vcf.gz

# Index the compressed VCF file
tabix -p vcf results_numenius/final_variants.clean.max_missing_0.8.recode.vcf.gz

# Second, filter with --max-missing 0.9 (keep variants with 10% or less missing data - stricter)
# Using --stdout to pipe directly to bgzip for compression
vcftools --gzvcf results_numenius/final_variants.clean.vcf.gz --keep results_numenius/Selected_SRR.tsv --max-missing 0.9 --recode --recode-INFO-all --stdout | bgzip -c > results_numenius/final_variants.clean.max_missing_0.9.recode.vcf.gz

# Index the compressed VCF file
tabix -p vcf results_numenius/final_variants.clean.max_missing_0.9.recode.vcf.gz

# Verify output files were created
ls -lh results_numenius/final_variants.clean.max_missing_0.*.recode.vcf*

# New filter MAT 
vcftools --gzvcf results_numenius/final_variants.clean.max_missing_0.9.recode.vcf.gz \
  --maf 0.05 \
  --recode --stdout | bgzip -c > results_numenius/output.maf05.vcf.gz

  # Index it
  tabix -p vcf results_numenius/output.maf05.vcf.gz

  ##  # Analyse and estimate inbreeding rate with PLINK ##
# Install plink
  cd ~
mkdir plink
cd plink
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip #linux version

unzip plink_linux_x86_64_20231211.zip
chmod +x plink

# Test to see if plink is loaded
./plink --help
export PATH=$PATH:~/plink
plink --help

# Set working directory
# Important to change to your working directory
cd /home/ba-student2/Personal_project_bioinformatic_applications/results_numenius

# PLINK Pipeline

plink --vcf output.maf05.vcf.gz \
--allow-extra-chr \
--make-bed \
--out data

plink --bfile data \
--allow-extra-chr \
--indep-pairwise 50 5 0.2 \
--out prune

plink --bfile data \
--allow-extra-chr \
--extract prune.prune.in \
--make-bed \
--out data_pruned

# Output PCA

plink --bfile data_pruned \
--allow-extra-chr \
--pca \
--out PCA

## Now to have the PCA figure, we use the results of PLINK (PCA.eigenvec) in RStudio. See the code called "PCA_figure.txt".
