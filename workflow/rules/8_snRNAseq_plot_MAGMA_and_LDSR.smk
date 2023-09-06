rule plot_MAGMA_and_LDSR_data:

    input:  magma: "../results/magma/snRNAseq_{GWAS}.herring.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out",
            ldsr: "../results/LDSR_part_herit/baseline_v1.2/herring/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    output: "../results/figures/{GWAS}_magma_ldsr_herring_lvl{LEVEL}_plot.jpeg"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: study_id = "herring",
            magma_dir = "../results/magma/",
            ldsr_dir = "../results/LDSR_part_herit/baseline_v1.2/herring/",
            fig_dir = "../results/figures/",
            level = "{LEVEL}"
    log:    "../results/00LOG/plot_magma_and_ldsr/snRNAseq_plot_magma_and_ldsr_{GWAS}_herring_lvl{LEVEL}.log"
    script:
            "../scripts/snRNAseq_plot_MAGMA_and_LDSR_new.R"

rule plot_MAGMA_and_LDSR_data_dwnSmpl:

    input:  magma: "../results/magma/snRNAseq_{GWAS}.herring_dwnSmpl.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out",
            ldsr: "../results/LDSR_part_herit/baseline_v1.2/herring_dwnSmpl/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    output: "../results/figures/{GWAS}_magma_ldsr_herring_dwnSmpl_lvl{LEVEL}_plot.jpeg"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: study_id = "herring_dwnSmpl",
            magma_dir = "../results/magma/",
            ldsr_dir = "../results/LDSR_part_herit/baseline_v1.2/herring_dwnSmpl/",
            fig_dir = "../results/figures/",
       	    level = "{LEVEL}"
    log:    "../results/00LOG/plot_magma_and_ldsr/snRNAseq_plot_magma_and_ldsr_{GWAS}_herring_dwnSmpl_lvl{LEVEL}.log"
    script:
            "../scripts/snRNAseq_plot_MAGMA_and_LDSR_new.R"

rule plot_MAGMA_and_LDSR_data_top2000:

    input:  magma: "../results/magma/snRNAseq_{GWAS}.herring_top2000.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out",
            ldsr: "../results/LDSR_part_herit/baseline_v1.2/herring/top2000_genes/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    output: "../results/figures/{GWAS}_magma_ldsr_herring_top2000_lvl{LEVEL}_plot.jpeg"
    resources: slurm_extra = "--use-singularity"
    singularity: "../resources/containers/snrna-seq_herring_complete_latest.sif"
    params: study_id = "herring_top2000",
            magma_dir = "../results/magma/",
            ldsr_dir = "../results/LDSR_part_herit/baseline_v1.2/herring/top2000_genes/",
            fig_dir = "../results/figures/",
       	    level = "{LEVEL}"
    log:    "../results/00LOG/plot_magma_and_ldsr/snRNAseq_plot_magma_and_ldsr_{GWAS}_herring_top2000_lvl{LEVEL}.log"
    script:
            "../scripts/snRNAseq_plot_MAGMA_and_LDSR_new.R"
