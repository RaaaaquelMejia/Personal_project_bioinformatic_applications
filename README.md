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

## Bibliography

Croll, D. (2026). genomepanel_nf - a highly efficient Nextflow pipeline for reference genome variant calling of large genome panels (v1.0.5). Zenodo. https://doi.org/10.5281/zenodo.19392838

Tan, H.Z., Ng, E.Y.X., Tang, Q. et al. Population genomics of two congeneric Palaearctic shorebirds reveals differential impacts of Quaternary climate oscillations across habitats types. Sci Rep 9, 18172 (2019). https://doi.org/10.1038/s41598-019-54715-9
