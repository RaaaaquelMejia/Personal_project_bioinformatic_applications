# Personal project in Bioinformatics Application - Raquel Mejia
## Estimation of inbreeding coefficients of the endangered bird species *Numenius arquata*

The Eurasian Curlew (*Numenius arquata*, Scolopacidae) is a widespread migratory species of wader found throughout temperate Europe, Asia and Africa.
In winter, the species is found in interdidal areas but in summer, for the breeding period, it migrates in open areas such as wet meadows or moorlands where it nests.

The loss of wetlands and the growing pression from agriculture are having an important impact on the species populations, which are declining. Since 2017, the species is assessed by IUCN as Near Threatened globaly (UICN, 2017), with an estimated global population of 610,000 to 830,000 mature individuals (Wetlands International, 2025). In Europe, the species is also assessed as Near Threatened (UICN, 2020), with the population estimated at 470,000 mature individuals.

The declin in population size is likely to have an impact on the genetic diversity of the species and the level of inbreeding among individuals. Therefore, the main goal of this project is to estimate the inbreeding coefficient of the Eurasian curlew using genetic data collected from 15 individuals across the World. 

# Methods

## Uploading genetic data

We used the available data of Tan et al. (2019), published with the accession number **PRJNA562783** in NCBI. The project gives access to 54 different SRA of *Numenius sp.* individuals. 
The complete list of SRA loaded is available here : https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=MCID_69f87e26a2960a7457636ada&o=acc_s%3Aa. We copied the run accession numbers into a text file called SRA.txt.

The reference genome used for the analysis (*Table 1*) was obtained from a *Numenius arquata* individual. It was dowloaded from NCBI web site using the accession number **PRJEB75987**.

<ins>Table 1 : Reference genome used for analysis </ins>
| Scientific name  | Assembly name | Assembly Accession | Source | Annotation | Level | Gene count |Contig N50|Size|BioProject|
| :---         |          ---: | :---         |     :---:      |          ---: |:---|:---|:---|:---|:---|
| *Numenius arquata*   | bNumArq3.hap1.1    | GCA_964106895.1  | GenBank    |     |Chromosome| |2605715|1348859203|PRJEB75987|
| *Numenius arquata*     | bNumArq3.hap1.1      | GCF_964106895.1  | RefSeq     | GCF_964106895.1-RS_2025_05    |Chromosome|17452|2605715|1348859203|PRJEB75987|

## Quality control and reference genome mapping

We used the genomepanel_nf Nextflow pipeline (Croll, 2026) to process the data. 

First of all, we installed genomepanel_nf directly on Visual Code Studio thanks to LEGCompute (Croll lab access needed). If you do not have access to Croll lab, you can load the pipeline here : https://crolllab.github.io/genomepanel_nf/getting-started/.

```
module load genomepanel_nf # Allows to install the pipeline
```
Next, we loaded the reference genome data and simplified it into a smaller FASTA substet containing chromosomals contings.

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/106/895/GCF_964106895.1_bNumArq3.hap1.1/GCF_964106895.1_bNumArq3.hap1.1_genomic.fna.gz
gunzip GCF_964106895.1_bNumArq3.hap1.1_genomic.fna.gz

# subsetting reference genome to chromosomal contigs (starting with ">NC_...")
module load SAMtools

samtools faidx GCF_964106895.1_bNumArq3.hap1.1_genomic.fna
grep "^>NC_" GCF_964106895.1_bNumArq3.hap1.1_genomic.fna | awk '{print substr($1,2)}' | xargs samtools faidx GCF_964106895.1_bNumArq3.hap1.1_genomic.fna > GCF_964106895.1_bNumArq3.hap1.1_genomic.subset.fna
```
After that, we were able to run the pipeline.

```
# run pipeline
screen -S gp

  nextflow run $GENOMEPANEL_HOME/pipeline/main.nf \
     --reference /home/ba-student2/Personal_project_bioinformatic_applications/GCF_964106895.1_bNumArq3.hap1.1_genomic.subset.fna \
     --reference_segments 20000000 \
     --SRA_index /home/ba-student2/Personal_project_bioinformatic_applications/SRA.txt \
     --ploidy 2 \
     --outdir results_numenius -profile slurm

screen -r gp
```
The results of the analysis were gived in a HTML file available as "pipeline_report_Numenius_arquata.html" in the main repository. 

## Data filtering to keep only the good quality data with VCFtools

VCFtools is a program package designed for working with VCF files. This toolset can be used to perform different operations on VCF files such as filtering out specific variants, create intersections and subsets of variants or compare files.

- Based on the results of the pipeline, we filtered to keep only the good data (with high concordance with the reference genome ?)
- We selected only the data with > 80 or 90 concordance and create a new file containing only this data

```
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
```
## Analyse and estimate inbreeding rate with PLINK

PLINK is an open-source whole genome association analysis toolset, designed to perform a range of basic, large-scale analyses in a computationally efficient manner.
We used PLINK here to analyse the data and generate a PCA output of the proximity of individuals based on their genetic data.

```
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
```

## Estimating the inbreeding rate with PCA visualisation

The PLINK pipeline created a nex file with the data needed to visualise a PCA of the individuals genetic data based on their species and location.
This data is contained in the file called "PCA_figure.txt". 
To visualise the PCA, we used ggplot2 on RStudio. The results are available in the main folder.

```
library(ggplot2)

pca <- read.table("PCA.eigenvec.txt", header=FALSE)

colnames(pca)[1:2] <- c("FID", "IID")

colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_minimal() +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA - Numenius arquata")


meta <- read.table("meta_2.txt", header = TRUE)

pca <- merge(pca, meta, by = "IID")

ggplot(pca, aes(x = PC1, y = PC2, color = Pop)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA - Numenius arquata and Numenius phaeopus",
    x = "PC1",
    y = "PC2",
    color = "Population"
  )

ggplot(pca, aes(x = PC1, y = PC2, color = Pop)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA - Numenius arquata",
    x = "PC1",
    y = "PC2",
    color = "Population"
  )

ggplot(pca, aes(PC1, PC2, color = Pop)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()
```

## Bibliography

Croll, D. (2026). genomepanel_nf - a highly efficient Nextflow pipeline for reference genome variant calling of large genome panels (v1.0.5). Zenodo. https://doi.org/10.5281/zenodo.19392838

Tan, H.Z., Ng, E.Y.X., Tang, Q. et al. Population genomics of two congeneric Palaearctic shorebirds reveals differential impacts of Quaternary climate oscillations across habitats types. Sci Rep 9, 18172 (2019). https://doi.org/10.1038/s41598-019-54715-9

The Variant Call Format and VCFtools, Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert Handsaker, Gerton Lunter, Gabor Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin and 1000 Genomes Project Analysis Group, Bioinformatics, 2011

Weeks JP (2010). “plink: An R Package for Linking Mixed-Format Tests Using IRT-Based Methods.” Journal of Statistical Software, 35(12), 1–33. http://www.jstatsoft.org/v35/i12/.
