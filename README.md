# README: Bulk RNA-Seq Analytics

## Overview

This repository contains a complete workflow for the analysis of bulk RNA-Seq data. It includes scripts and documentation for processing, analyzing, and visualizing high-throughput sequencing (HTS) data. The project is focused on comparing the performance of differential expression analysis pipelines (`edgeR` and `limma-voom`) using RNA-Seq data from a publicly available dataset.

## Key Features

- Programmatic data acquisition from the NCBI GEO database.
- Preprocessing steps including quality control, trimming, and alignment.
- Comprehensive differential expression analysis with `edgeR` and `limma-voom`.
- Visualization of results including MA plots, MDS plots, and word clouds.
- Analysis of under-sampled and full-count datasets for comparative evaluation.

## Dataset

The RNA-Seq data analyzed in this project is from:

> **Fu, N.Y. et al. (2015).** EGF-mediated induction of *Mcl-1* at the switch to lactation is essential for alveolar cell survival. [Nature Cell Biology](https://www.nature.com/articles/ncb3117).  
> GEO Accession: [GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450)

## Prerequisites

- **Programming Languages:** R (≥ 4.0), UNIX shell scripting
- **R Packages:**
  - `edgeR`
  - `limma`
  - `Rsubread`
  - `org.Mm.eg.db`
  - `DT`
  - `wordcloud2`
- **Tools:**
  - [SRA Toolkit](https://github.com/ncbi/sra-tools)
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [MultiQC](https://multiqc.info/)
  - [fastp](https://github.com/OpenGene/fastp)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/username/bulk-rna-seq.git
   cd bulk-rna-seq
   ```
2. Install required R packages:
   ```r
   install.packages("BiocManager")
   BiocManager::install(c("edgeR", "limma", "Rsubread", "org.Mm.eg.db", "DT", "wordcloud2"))
   ```
3. Install additional tools (`SRA Toolkit`, `fastp`, `FastQC`, `MultiQC`) using your system's package manager.
4. Obtain reference genome and its associated annotation data (see below).

## Outputs

- **Differentially Expressed Genes (DEGs)**: Lists of DEGs with $p$-values and fold changes.
- **Quality Control Reports**: Aggregated results from `FastQC` and `MultiQC`.
- **Plots**:
  - MA plots and MDS plots
  - Word clouds of DEGs
  - Interactive Tables: Sortable and searchable tables of DEGs.

## Project Goals

This project aims to:
- Develop a reproducible workflow for bulk RNA-Seq data analysis.
- Compare the performance of `edgeR` and `limma-voom` pipelines.
- Highlight the challenges and considerations in RNA-Seq data analysis.

## Reference Genome and Annotation Data

In your project directory, create a new folder `Mus_musculus.GRCm39`

Navigate over to [Ensembl genome browser 113] (<http://useast.ensembl.org/index.html>) and Find the Page for the Mouse (GRCm39) Reference Assembly. Although there are many "images" or "version" of this assembly that include or exclude various products of the assembly process, a good starting point is often the "primary assembly" that's been "masked" for repetitive elements and low-complexity regions.

Click on the link "Download DNA sequence (FASTA)" to reach an FTP-page with genome data. According to the `README` file at the bottom of the file-list:

```         
<sequence type>:
  * 'dna' - unmasked genomic DNA sequences.
  * 'dna_rm' - masked genomic DNA.  Interspersed repeats and low
     complexity regions are detected with the RepeatMasker tool and masked
     by replacing repeats with 'N's.
  * 'dna_sm' - soft-masked genomic DNA. All repeats and low complexity regions
    have been replaced with lowercased versions of their nucleic base
```

We will work with the masked primary assembly for this assignment, so please download the file `Mus_musculus.GRCm39.dna_rm.primary_assembly.fa.gz`.

Back out to the main page for "Mouse (GRCm39)", in the Gene Annotation panel on the top-right, get the associated annotation data paired to the reference genome assembly in Release 113 in GTF format by clicking "Download GTF" to obtain the file `Mus_musculus.GRCm39.113.gtf.gz`

Inside your project directory for this lab, move your Ensembl genome data files to the subdirectory `Mus_musculus.GRCm39`.

After downloading binaries or compiling, you will only need to install `subread-buildindex`, `subread-align`, and `featureCounts`, for example by copying these executable binaries to `/usr/bin/local` and possibly also re-starting your shell or terminal.

From the project directory (parent directory to where the genome data now is) do:

```bash         
subread-buildindex -o Mus_musculus.GRCm39/Mus_musculus.GRCm39_subread Mus_musculus.GRCm39/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa.gz
```


