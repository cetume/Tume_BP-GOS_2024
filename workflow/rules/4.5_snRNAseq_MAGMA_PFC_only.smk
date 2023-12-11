rule magma_map_snps_to_genes_PFC_only:
    input:   snp_loc = "../resources/refs/g1000_eur.bim",
             gene_loc = "../results/R_objects/Ensembl_PFC_genes.hg19.extendedMHCremoved.gene.loc.txt"
    output:  "../results/magma/snRNAseq_PFC_only.magma.35UP_10DOWN.genes.annot"
    params:  "../results/magma/snRNAseq_PFC_only.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma annotation step to map SNPs to genes."
    log:     "../results/logs/magma/snRNAseq_PFC_only.annotate.snps2genes.35UP_10DOWN.log"
    run:

            shell("""

            module load magma/1.10
            magma --annotate window=35,10 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}

            """)

rule magma_gene_analysis_PFC_only:
    input:   gene_annot = "../results/magma/snRNAseq_PFC_only.magma.35UP_10DOWN.genes.annot",
             gwas = "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    output:  "../results/magma/snRNAseq_{GWAS}_PFC_only.magma.35UP_10DOWN.genes.raw",
             "../results/magma/snRNAseq_{GWAS}_PFC_only.magma.35UP_10DOWN.genes.out"
    params:  ref = "../resources/refs/g1000_eur",
             out = "../results/magma/snRNAseq_{GWAS}_PFC_only.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene analysis step for {wildcards.GWAS}"
    log:     "../results/logs/magma/snRNAseq_PFC_only.gene_analysis.{GWAS}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --bfile {params.ref} --pval {input.gwas} ncol='N' --gene-annot {input.gene_annot} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis_PFC_only:
    input:   genes = "../results/magma/snRNAseq_{GWAS}_PFC_only.magma.35UP_10DOWN.genes.raw",
             data = "../results/gene_lists/herring/MAGMA/herring_lvl{LEVEL}.txt"
    output:  "../results/magma/snRNAseq_{GWAS}_PFC_only.herring.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_{GWAS}_PFC_only.herring.lvl{LEVEL}.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for {wildcards.GWAS}, cluster level {wildcards.LEVEL}"
    log:     "../results/logs/magma/snRNAseq_PFC_only.gene_set_analysis.{GWAS}.herring.35UP_10DOWN.lvl{LEVEL}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """


