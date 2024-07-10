# Genetic Implication of Specific Glutamatergic Neurons of the Prefrontal Cortex in the Pathophysiology of Schizophrenia

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). The workflow follows the the snakemake distribution and reproducibility recommendations.

A snakemake pipeline to process snRNAseq data, from downloading the data to testing gene set enrichment.

The snRNAseq data for this study was taken from [Herring et al. (2022)](https://doi.org/10.1016/j.cell.2022.09.039):

* [Matrix and metadata](https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_full-counts-and-downsampled-CPM.h5ad)

#### Scripts:
1. [snRNAseq_prep_herring_data.py](workflow/scripts/snRNAseq_prep_herring_data.py) - Extract cell and gene metadata from snRNAseq dataset
2. [snRNAseq_prep_herring_data.R](workflow/scripts/snRNAseq_prep_herring_data.R) - Prepare Herring data for analysis 
3. [snRNAseq_rm_MHC_from_ref.R](workflow/scripts/snRNAseq_rm_MHC_from_ref.R) - Prepare genome reference for MAGMA and LDSR
4. [snRNAseq_prep_enrich_files.R](workflow/scripts/snRNAseq_prep_enrich_files.R) - Generate specificity scores and top 10% genesets for enrichment analyses
5. [snRNAseq_prep_enrich_files_only.R](workflow/scripts/snRNAseq_prep_enrich_files_only.R) - Generate top 10% genesets for enrichment analyses when ctd object is already available
<!-- 6. snRNAseq_plot_MAGMA_and_LDSR.R -->

#### Singularity Container:
The singularity container used to run the R scripts in this pipeline was generated in [Sylabs](https://cloud.sylabs.io/), using the following definition file:

```
Bootstrap: docker
From: bioconductor/bioconductor_docker:devel

%labels
    Version v0.0.1

%help
    This is a container to run the analysis for the herring 2022 single nuclei RNA seq data

%post
    R --no-echo -e 'install.packages("scTenifoldNet")'
    R --no-echo -e 'install.packages("dplyr")'
    R --no-echo -e 'install.packages("readr")'
    R --no-echo -e 'BiocManager::install("biomaRt")'
    R --no-echo -e 'BiocManager::install("EWCE")'
    R --no-echo -e 'install.packages("ggplot2")'
    R --no-echo -e 'install.packages("Seurat")'
    R --no-echo -e 'install.packages("Matrix")'
    R --no-echo -e 'install.packages("ggsignif")'
    R --no-echo -e 'install.packages("cowplot")'
    R --no-echo -e 'install.packages("reshape2")'
    R --no-echo -e 'install.packages("ggdendro")'
    R --no-echo -e 'install.packages("data.table")'
    R --no-echo -e 'install.packages("readxl")'
    R --no-echo -e 'install.packages("scCustomize")'
    R --no-echo -e 'BiocManager::install("SeuratWrappers")'
    R --no-echo -e 'BiocManager::install("scuttle")'
```

#### Running the Snakemake Pipeline:
The Snakemake file includes rules from 11 files. To run the entire pipeline, rules from files 1-3 must be run first before running rules from files 4-9, by blanking out the relevant lines in the Snakemake file using #.

By default, this pipeline will analyse 84 cell sub-cluster populations. 

#### Copyright and Licence Information:
See the [LICENCE file](LICENCE.md)
