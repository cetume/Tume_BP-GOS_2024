rule plot_MAGMA_and_LDSR_data:
    input:  magma: "../results/magma/snRNAseq_{GWAS}.{project}.lvl_{LEVEL}.magma.35UP_10DOWN.gsa.out",
            ldsr: "../results/LDSR_part_herit/baseline_v1.2/{project}/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    output: "../results/figures/{GWAS}_magma_ldsr_{project}_lvl_{LEVEL}_plot.jpeg"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: magma_dir = "../results/magma/",
            ldsr_dir = "../results/LDSR_part_herit/baseline_v1.2/{project}/",
            fig_dir = "../results/figures/",
            level = "{LEVEL}"
    log:    "../results/00LOG/plot_magma_and_ldsr/snRNAseq_plot_magma_and_ldsr_{project}_lvl{LEVEL}.log"
    script:
            "../scripts/snRNAseq_plot_MAGMA_and_LDSR.R"
