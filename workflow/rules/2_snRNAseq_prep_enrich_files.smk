localrules: rm_MHC_from_ref

## Remove MHC region and XY genes from reference file
## Outputs two files:
##  - the full ref file with MHC and XY genes removed
##  - R object for next rule of above with fewer cols

rule rm_MHC_from_ref:
    #Creates the reference files for gene locations - all genes, protein-coding genes only and non-coding genes only
    output: gene_coord = "../results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds",
            mhc_genes = "../results/R_objects/mhc_genes.rds",
            xy_genes = "../results/R_objects/xy_genes.rds",
            protein_gene_coord = "../results/R_objects/Ensembl_hg19_protein_gene_coords_noMHC.rds",
            nonprotein_gene_coord = "../results/R_objects/Ensembl_hg19_nonprotein_gene_coords_noMHC.rds"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: ref_noMHC = "../results/R_objects/Ensembl.hg19.MHCremoved.gene.loc.txt",
            protein_ref_noMHC = "../results/R_objects/Ensembl.hg19.MHCremoved.protein.gene.loc.txt",
            nonprotein_ref_noMHC = "../results/R_objects/Ensembl.hg19.MHCremoved.nonprotein.gene.loc.txt"
    threads: 1
    log:    "../results/00LOG/prep_enrich_files/snRNAseq_rm_MHC_from_ref.log"
    script:
            "../scripts/snRNAseq_rm_MHC_from_ref.R"

rule prep_enrichment_files:
    #Creates the ctd object and the MAGMA and SLDSR gene lists
    input:  seurat_obj = "../results/R_objects/seurat_{project}.rds",
            gene_coord = "../results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds",
            protein_gene_coord = "../results/R_objects/Ensembl_hg19_protein_gene_coords_noMHC.rds",
            nonprotein_gene_coord = "../results/R_objects/Ensembl_hg19_nonprotein_gene_coords_noMHC.rds",
            mhc_genes = "../results/R_objects/mhc_genes.rds",
            xy_genes = "../results/R_objects/xy_genes.rds"
    output: "../results/gene_lists/prep_enrichment_file_{project}.done"
    resources: tasks = 1, mem_mb = 240000, slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: ctd_outdir = "../results/ctd_objects/",
            outdir = "../results/gene_lists/",
            study_id = "{project}"
    threads: 10
    log:    "../results/00LOG/prep_enrich_files/snRNAseq_prep_enrichment_files_{project}.log"
    script:
            "../scripts/snRNAseq_prep_enrich_files_final.R"

rule prep_enrichment_files_only:
    #Can use to create the MAGMA and SLDSR gene lists if the ctd object already exists
    input:  ctd_object = "../results/ctd_objects/ctd_{project}.rda",
            gene_coord = "../results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds",
            protein_gene_coord = "../results/R_objects/Ensembl_hg19_protein_gene_coords_noMHC.rds",
            nonprotein_gene_coord = "../results/R_objects/Ensembl_hg19_nonprotein_gene_coords_noMHC.rds"
    output: "../results/gene_lists/prep_enrichment_file_only_{project}.done"
    resources: tasks = 1, mem_mb = 15000, slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: outdir = "../results/gene_lists/",
            study_id = "{project}"
    threads: 1
    log:    "../results/00LOG/prep_enrich_files/snRNAseq_prep_enrichment_files_{project}.log"
    script:
            "../scripts/snRNAseq_prep_enrich_files_only.R"
