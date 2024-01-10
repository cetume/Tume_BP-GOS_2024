rule magma_map_snps_to_genes:
    input:   snp_loc = "../resources/refs/g1000_eur.bim",
             gene_loc = "../results/R_objects/Ensembl.hg19.MHCremoved.gene.loc.txt"
    output:  "../results/magma/snRNAseq.magma.35UP_10DOWN.genes.annot"
    params:  "../results/magma/snRNAseq.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma annotation step to map SNPs to genes."
    log:     "../results/logs/magma/snRNAseq.annotate.snps2genes.35UP_10DOWN.log"
    run:

            shell("""

            module load magma/1.10
            magma --annotate window=35,10 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}

            """)

rule magma_gene_analysis:
    input:   gene_annot = "../results/magma/snRNAseq.magma.35UP_10DOWN.genes.annot",
             gwas = "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    output:  "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.genes.raw",
             "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.genes.out"
    params:  ref = "../resources/refs/g1000_eur",
             out = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene analysis step for {wildcards.GWAS}"
    log:     "../results/logs/magma/snRNAseq.gene_analysis.{GWAS}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --bfile {params.ref} --pval {input.gwas} ncol='N' --gene-annot {input.gene_annot} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis:
    input:   genes = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.genes.raw",
             data = "../results/gene_lists/herring/MAGMA/herring_lvl{LEVEL}.txt"
    output:  "../results/magma/snRNAseq_{GWAS}.herring.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_{GWAS}.herring.lvl{LEVEL}.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for {wildcards.GWAS}, cluster level {wildcards.LEVEL}"
    log:     "../results/logs/magma/snRNAseq.gene_set_analysis.{GWAS}.herring.35UP_10DOWN.lvl{LEVEL}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis_dwnSmpl:
    input:   genes = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.genes.raw",
             data = "../results/gene_lists/herring_dwnSmpl/MAGMA/herring_dwnSmpl_lvl{LEVEL}.txt"
    output:  "../results/magma/snRNAseq_{GWAS}.herring_dwnSmpl.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_{GWAS}.herring_dwnSmpl.lvl{LEVEL}.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for dwnSmpl {wildcards.GWAS}, cluster level {wildcards.LEVEL}"
    log:     "../results/logs/magma/snRNAseq.gene_set_analysis.{GWAS}.herring_dwnSmpl.35UP_10DOWN.lvl{LEVEL}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis_top2000:
    input:   genes = "../results/magma/snRNAseq_{GWAS}.magma.35UP_10DOWN.genes.raw",
             data = "../results/gene_lists/herring/MAGMA/herring_top2000_lvl{LEVEL}.txt"
    output:  "../results/magma/snRNAseq_{GWAS}.herring_top2000.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_{GWAS}.herring_top2000.lvl{LEVEL}.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for top2000 {wildcards.GWAS}, cluster level {wildcards.LEVEL}"
    log:     "../results/logs/magma/snRNAseq.gene_set_analysis.{GWAS}.top2000.35UP_10DOWN.lvl_{LEVEL}.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

rule magma_conditional:
    input:   gene_list = "../results/gene_lists/herring/MAGMA/herring_lvl{LEVEL}.txt",
             scz_magma = "../results/magma/snRNAseq_SCZ.magma.35UP_10DOWN.genes.raw" 
    output:  "../results/magma_conditional/snRNAseq.herring_all_sig_condition_{CONDITION}.lvl{LEVEL}.magma.35UP_10DOWN.gsa.out"
    params:  "../results/magma_conditional/snRNAseq.herring_all_sig_condition_{CONDITION}.lvl{LEVEL}.magma.35UP_10DOWN"
    message: "Running MAGMA on all significant cell types conditioning on {wildcards.CONDITION}, lvl{wildcards.LEVEL}"
    log:     "../results/logs/magma/magma_conditional/snRNAseq.magma.conditional.{CONDITION}.lvl{LEVEL}.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} --model condition={wildcards.CONDITION} --out {params} &> {log}

             """

rule magma_gene_set_analysis_GO:
    input:   genes = "../results/magma/snRNAseq_SCZ.magma.35UP_10DOWN.genes.raw",
             data = "../results/gene_lists/herring/MAGMA/GO_term_genes_for_magma.txt"
    output:  "../results/magma/snRNAseq_SCZ.GO_term_genes.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_SCZ.GO_term_genes.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for SCZ, GO term genes"
    log:     "../results/logs/magma/snRNAseq.gene_set_analysis.SCZ.herring_GO_term_genes.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

rule magma_GO_interactions:
    input:   genes = "../results/magma/snRNAseq_SCZ.magma.35UP_10DOWN.genes.raw",
             data = "../results/gene_lists/herring/MAGMA/GO_term_genes_for_cond_magma_test.txt",
             interaction_list = "../resources/sheets/interaction_list"
    output:  "../results/magma/snRNAseq_SCZ.GO_term_genes_inter.magma.35UP_10DOWN.gsa.out"
    params:  out = "../results/magma/snRNAseq_SCZ.GO_term_genes_inter.magma.35UP_10DOWN"
    resources: slurm_extra = "--use-conda"
    message: "Running magma gene set analysis step for SCZ, GO term genes interactions"
    log:     "../results/logs/magma/snRNAseq.gene_set_analysis.SCZ.herring_GO_term_genes_inter.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --model interaction={input.interaction_list} interaction-ss-size=25,0.03 --out {params.out} &> {log}

             """
