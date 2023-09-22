rule magma_map_snps_to_protein_genes:
    input:   snp_loc = "../resources/refs/g1000_eur.bim",
             gene_loc = "../results/R_objects/Ensembl.hg19.MHCremoved.protein.gene.loc.txt"
    output:  "../results/magma/snRNAseq.magma.35UP_10DOWN.protein.genes.annot"
    params:  "../results/magma/snRNAseq.magma.35UP_10DOWN.protein"
    resources: slurm_extra = "--use-conda"
    message: "Running magma annotation step to map SNPs to protein-coding genes."
    log:     "../results/logs/magma/snRNAseq.annotate.snps2proteingenes.35UP_10DOWN.log"
    run:

            shell("""

            module load magma/1.10
            magma --annotate window=35,10 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}

            """)

rule magma_gene_analysis:
    input:   gene_annot = "../results/magma/snRNAseq.magma.35UP_10DOWN.protein.genes.annot",
             gwas = "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    output:  "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.protein.genes.raw",
             "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.protein.genes.out"
    params:  ref = "../resources/refs/g1000_eur",
             out = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.protein"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene analysis step for {wildcards.GWAS}"
    log:     "../results/logs/magma/snRNAseq.protein_gene_analysis.{GWAS}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --bfile {params.ref} --pval {input.gwas} ncol='N' --gene-annot {input.gene_annot} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis:
    input:   genes = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.protein.genes.raw",
             data = "../results/gene_lists/herring/MAGMA/herring_protein_coding_lvl{LEVEL}.txt"
    output:  "../results/magma/snRNAseq_{GWAS}.herring_protein_coding.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_{GWAS}.herring_protein_coding.lvl{LEVEL}.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for {wildcards.GWAS}, cluster level {wildcards.LEVEL}, protein-coding genes"
    log:     "../results/logs/magma/snRNAseq.protein_gene_set_analysis.{GWAS}.herring.35UP_10DOWN.lvl{LEVEL}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """
