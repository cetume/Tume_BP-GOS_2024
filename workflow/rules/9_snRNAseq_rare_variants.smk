rule rare_variants:
    input:  ctd_object = "../results/ctd_objects/ctd_{project}.rda",
            schema = "../resources/rare_variants/meta_results_2023_05_30_15_17_13.csv",
            schema_genes = "../resources/rare_variants/schema_gene_conversion.txt"
    output: "../results/rare_variants/rare_variants_analysis_{project}.done",
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: study_id = "{project}",
            rare_dir = "../results/rare_variants/",
            boxplots_dir = "../results/rare_variants/boxplots_temp/{project}/",
            fig_dir = "../results/figures/",
            ctd_dir = "../results/ctd_objects/"
    log:    "../results/00LOG/rare_variants/snRNAseq_rare_variants_{project}.log"
    script:
            "../scripts/snRNAseq_rare_variants_new.R"
