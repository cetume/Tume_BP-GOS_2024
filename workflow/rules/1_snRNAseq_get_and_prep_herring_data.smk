localrules: download_herring_data

## Download Herring et al 2022 data
## Data from https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_full-counts-and-downsampled-CPM.h5ad
## Generate count and meta.data files for R

rule download_herring_data:
    output:  "../resources/raw_data/herring_2022/RNA-all_full-counts-and-downsampled-CPM.h5ad"
    params:  outdir = "../resources/raw_data/herring_2022/"
    message: "Downloading herring data"
    log:     "../results/logs/herring/snRNAseq.download.herring.log"
    shell:
             """
             wget https://storage.googleapis.com/neuro-dev/Processed_data/RNA-all_full-counts-and-downsampled-CPM.h5ad -P {params.outdir}
             """

rule prep_herring_data_for_R:
    input:   "../resources/raw_data/herring_2022/RNA-all_full-counts-and-downsampled-CPM.h5ad"
    output:  counts = "../resources/raw_data/herring_2022/data_for_R/herring_counts.mtx",
             cellMeta = "../resources/raw_data/herring_2022/data_for_R/herring_counts_cellMeta.csv",
             geneMeta = "../resources/raw_data/herring_2022/data_for_R/herring_counts_geneMeta.csv"
    resources: tasks = 1, mem_mb = 15000, slurm_extra = "--use-conda"
    message: "Preparing herring data for R"
    log:     "../results/logs/herring/snRNAseq.prep.herring.log"
    script:  
             "../scripts/snRNAseq_prep_herring_data.py"

rule prep_herring_data_and_dwnSmpl_for_ctd:
    input:  counts = "../resources/raw_data/herring_2022/data_for_R/herring_counts.mtx",
            cellMeta = "../resources/raw_data/herring_2022/data_for_R/herring_counts_cellMeta.csv",
            geneMeta = "../resources/raw_data/herring_2022/data_for_R/herring_counts_geneMeta.csv"
    output: herring = "../results/R_objects/seurat_herring.rds", 
            herring_dwnSmpl_2 = "../results/R_objects/seurat_herring_dwnSmpl_lvl2.rds"
    resources: tasks = 1, mem_mb = 25000, slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: r_dir = "../results/R_objects/"
    log:    "../results/00LOG/prep_enrich_files/snRNAseq_prep_herring_data.log"
    script:
            "../scripts/snRNAseq_prep_herring_data_and_dwnSmpl_final.R"
