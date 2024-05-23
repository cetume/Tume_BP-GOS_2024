rule magma_map_snps_to_genes_nonprotein_coding:
    input:   snp_loc = "../resources/refs/g1000_eur.bim",
             gene_loc = "../results/R_objects/Ensembl.hg19.MHCremoved.nonprotein.gene.loc.txt"
    output:  "../results/magma/snRNAseq.magma.35UP_10DOWN.nonprotein.genes.annot"
    params:  "../results/magma/snRNAseq.magma.35UP_10DOWN.nonprotein"
    resources: slurm_extra = "--use-conda"
    message: "Running magma annotation step to map SNPs to nonprotein-coding genes."
    log:     "../results/logs/magma/snRNAseq.annotate.snps2nonproteingenes.35UP_10DOWN.log"
    run:

            shell("""

            module load magma/1.10
            magma --annotate window=35,10 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}

            """)

rule magma_gene_analysis_nonprotein_coding:
    input:   gene_annot = "../results/magma/snRNAseq.magma.35UP_10DOWN.nonprotein.genes.annot",
             gwas = "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    output:  "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.nonprotein.genes.raw",
             "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.nonprotein.genes.out"
    params:  ref = "../resources/refs/g1000_eur",
             out = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.nonprotein"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene analysis step for {wildcards.GWAS}"
    log:     "../results/logs/magma/snRNAseq.nonprotein_gene_analysis.{GWAS}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --bfile {params.ref} --pval {input.gwas} ncol='N' --gene-annot {input.gene_annot} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis_nonprotein_coding:
    input:   genes = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.nonprotein.genes.raw",
             data = "../results/gene_lists/herring/MAGMA/herring_nonprotein_coding_lvl{LEVEL}.txt"
    output:  "../results/magma/snRNAseq_{GWAS}.herring_nonprotein_coding.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_{GWAS}.herring_nonprotein_coding.lvl{LEVEL}.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for {wildcards.GWAS}, cluster level {wildcards.LEVEL}, nonprotein-coding genes"
    log:     "../results/logs/magma/snRNAseq.nonprotein_gene_set_analysis.{GWAS}.herring.35UP_10DOWN.lvl{LEVEL}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """
