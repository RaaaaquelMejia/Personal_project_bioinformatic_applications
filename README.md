# Personal project in Bioinformatics Application - Raquel Mejia
## Estimation of inbreeding coefficients of the endangered bird species *Numenius arquata*

The Eurasian Curlew (*Numenius arquata*, Scolopacidae) is a widespread migratory species of wader found throughout temperate Europe, Asia and Africa.
In winter, the species is found in interdidal areas but in summer, for the breeding period, it migrates in open areas such as wet meadows or moorlands where it nests.

The loss of wetlands and the growing pression from agriculture are having an important impact on the species populations, which are declining. Since 2017, the species is assessed by IUCN as Near Threatened globaly (UICN, 2017), with an estimated global population of 610,000 to 830,000 mature individuals (Wetlands International, 2025). In Europe, the species is also assessed as Near Threatened (UICN, 2020), with the population estimated at 470,000 mature individuals.

The declin in population size is likely to have an impact on the genetic diversity of the species and the level of inbreeding among individuals. Therefore, the main goal of this project is to estimate the inbreeding coefficient of the Eurasian curlew using genetic data collected from 15 individuals across the World. 

## Genetic data used in the project

We used the available data of Tan et al. (2019), published with the accession number **PRJNA562783** in NCBI. The project gives access to 54 different SRA of *Numenius sp.* individuals. The complete list of SRA is available here : https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=MCID_69f87e26a2960a7457636ada&o=acc_s%3Aa

**Source :** Tan, H.Z., Ng, E.Y.X., Tang, Q. et al. Population genomics of two congeneric Palaearctic shorebirds reveals differential impacts of Quaternary climate oscillations across habitats types. Sci Rep 9, 18172 (2019). https://doi.org/10.1038/s41598-019-54715-9

For the reference genome, we were able to use a genome from a *Numenius arquata* individual (*Table 1*).

<ins>Table 1 : Reference genome used for analysis </ins>
| Scientific name  | Assembly name | Assembly Accession | Source | Annotation | Level | Gene count |Contig N50|Size|BioProject|
| :---         |          ---: | :---         |     :---:      |          ---: |:---|:---|:---|:---|:---|
| *Numenius arquata*   | bNumArq3.hap1.1    | GCA_964106895.1  | GenBank    |     |Chromosome| |2605715|1348859203|PRJEB75987|
| *Numenius arquata*     | bNumArq3.hap1.1      | GCF_964106895.1  | RefSeq     | GCF_964106895.1-RS_2025_05    |Chromosome|17452|2605715|1348859203|PRJEB75987|

## Methods

To process the data, we used the genomepanel_nf Nextflow pipeline (2026).
Source : Croll, D. (2026). genomepanel_nf - a highly efficient Nextflow pipeline for reference genome variant calling of large genome panels (v1.0.5). Zenodo. https://doi.org/10.5281/zenodo.19392838

```
module load genomepanel_nf # Allows to install the pipeline

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/106/895/GCF_964106895.1_bNumArq3.hap1.1/GCF_964106895.1_bNumArq3.hap1.1_genomic.fna.gz
gunzip GCF_964106895.1_bNumArq3.hap1.1_genomic.fna.gz

# subsetting reference genome to chromosomal contigs (starting with ">NC_...")
module load SAMtools

samtools faidx GCF_964106895.1_bNumArq3.hap1.1_genomic.fna
grep "^>NC_" GCF_964106895.1_bNumArq3.hap1.1_genomic.fna | awk '{print substr($1,2)}' | xargs samtools faidx GCF_964106895.1_bNumArq3.hap1.1_genomic.fna > GCF_964106895.1_bNumArq3.hap1.1_genomic.subset.fna


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
