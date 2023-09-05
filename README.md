# Investigating Cells Mediating Schizophrenia Genetic Risk using Single-Nuclei RNA Sequencing Data from Herring et al. (2022)

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). The paper is here. The workflow follows the the snakemake distribution and reproducibility recommendations.

A snakemake pipeline to process snRNAseq data from downloading the snRNAseq data to gene set enrichment analyses.

The snRNAseq data for this study was taken from [Herring et al. (2022)](https://doi.org/10.1016/j.cell.2022.09.039):

* [Matrix and metadata](https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_full-counts-and-downsampled-CPM.h5ad)

#### Scripts:
1. snRNAseq_prep_herring_data.py - Extract cell and gene metadata from dataset
2. snRNAseq_prep_herring_data_and_dwnSmpl.R - Prepare Herring data for analysis 
3. snRNAseq_rm_MHC_from_ref_new.R - Prepare genome reference for MAGMA and LDSR (identify and remove genes in MHC region)
4. snRNAseq_prep_enrich_files_methodB.R - Generate specificity scores and top 10% genesets for enrichment analyses
5. snRNAseq_get_and_munge_GWAS_for_MAGMA.sh - Download and prepare GWAS summary statistics for [schizophrenia](https://figshare.com/ndownloader/files/28169757) ([Trubetskoy et al., 2022](https://doi.org/10.1038/s41586-022-04434-5)) and [height](https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz) ([Yengo et al., 2018](https://doi.org/10.1093/hmg/ddy271))
<!-- 6. snRNAseq_plot_MAGMA_and_LDSR.R -->
<!-- 7. snRNAseq_plot_rare_variants.R -->

#### Running Snakemake Pipeline:
The Snakemake file includes rules from 9 files which encompass the entire analysis pipeline. To run the entire pipeline, rules from files 1-3 must be run first before running rules from files 4-9, by blanking out the relevant lines in the Snakemake file using #.

By default, this pipeline will analysis the dataset using two annotation levels:
* Level 1 = 19 major cell clusters
* Level 2 = 84 sub cell clusters

If you choose to run only annotation level 1 or level 2, the cell types not included must be blanked out and the required celltypes.tsv (both, lvl1_only, lvl2_only) file for LDSR must be select in the config.yaml file.
