rule rare_variants:
    input:  ctd_object = "../results/ctd_objects/ctd_{project}.rda",
            schema = "../resources/rare_variants/meta_results_2023_05_30_15_17_13.csv",
            schema_genes = "../resources/rare_variants/schema_gene_conversion.txt"
    output: "../results/rare_variants/rare_variants_analysis_{project}.done",
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: study_id = "{project}",
            rare_dir = "../results/rare_variants/",
            fig_dir = "../results/rare_variants/boxplots_temp/{project}/",
            ctd_dir = "../results/ctd_objects/"
    log:    "../results/00LOG/rare_variants/snRNAseq_rare_variants_{project}.log"
    script:
            "../scripts/snRNAseq_rare_variants.R"

rule plot_rare_variants:
    input:  "../results/rare_variants/wilcoxon_df_herring_lvl{LEVEL}.txt"
    output: "../results/figures/wilcoxon_herring_lvl{LEVEL}_plot.jpeg"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: fig_dir = "../results/figures/",
            level = "{LEVEL}"
    log:    "../results/00LOG/rare_variants/snRNAseq_plot_rare_variants_{LEVEL}.log"
    script:
            "../scripts/snRNAseq_plot_rare_variants_new.R"

rule plot_rare_variants_dwnSmpl:
    input:  "../results/rare_variants/wilcoxon_df_herring_dwnSmpl_lvl{LEVEL}.txt"
    output: "../results/figures/wilcoxon_herring_dwnSmpl_lvl{LEVEL}_plot.jpeg"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: fig_dir = "../results/figures/",
            level = "{LEVEL}"
    log:    "../results/00LOG/rare_variants/snRNAseq_plot_rare_variants_dwnSmpl_{LEVEL}.log"
    script:
            "../scripts/snRNAseq_plot_rare_variants_dwnSmpl_new.R"
