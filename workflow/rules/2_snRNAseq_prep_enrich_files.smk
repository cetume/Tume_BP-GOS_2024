localrules: rm_MHC_from_ref

## Remove MHC region genes from reference file
## Outputs two files:
##  - the full ref file with MHC genes removed
##  - R object for next rule of above with fewer cols

rule rm_MHC_from_ref:
    output: gene_coord = "../results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds",
            mhc_genes = "../results/R_objects/mhc_genes.rds"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: ref_noMHC = "../results/R_objects/Ensembl.hg19.MHCremoved.gene.loc.txt",
    threads: 1
    log:    "../results/00LOG/prep_enrich_files/snRNAseq_rm_MHC_from_ref.log"
    script:
            "../scripts/snRNAseq_rm_MHC_from_ref_new.R"

rule prep_enrichment_files:
    input:  seurat_obj = "../results/R_objects/seurat_{project}.rds",
            gene_coord = "../results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds",
            mhc_genes = "../results/R_objects/mhc_genes.rds"
    output: "../results/gene_lists/prep_enrichment_file_{project}.done"
    resources: tasks = 1, mem_mb = 240000, slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: ctd_outdir = "../results/ctd_objects/",
            outdir = "../results/gene_lists/",
            study_id = "{project}",
    threads: 10
    log:    "../results/00LOG/prep_enrich_files/snRNAseq_prep_enrichment_files_{project}.log"
    script:
            "../scripts/snRNAseq_prep_enrich_files_methodB.R"

