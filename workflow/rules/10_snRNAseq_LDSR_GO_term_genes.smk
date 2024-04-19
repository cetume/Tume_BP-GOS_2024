localrules: ldsr_stratified_summary_go_terms

rule ldsr_make_annot_go_terms:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   gene_set = "../results/gene_lists/herring/LDSR_GO_term_genes/{CELL_TYPE_GO_TERM}.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/LDSR_annotation_files/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{CHR}.annot.gz"
    resources: slurm_extra = "--use-conda"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snRNAseq: {wildcards.CELL_TYPE_GO_TERM}, Chr {wildcards.CHR}"
    log:     "../results/00LOG/LDSR/herring_GO_term_genes/make_annot.snRNAseq.{CELL_TYPE_GO_TERM}.Chr{CHR}.log"
    shell:
             """

             python ../resources/ldsr/make_annot.py \
             --bed-file {input.gene_set} \
             --windowsize 0 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}

             """

rule ldsr_ld_scores_go_terms:
    input:   annot = "../results/LDSR_annotation_files/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsr/reference_files/hapmap3_snps"
    output:  "../results/LDSR_annotation_files/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{CHR}.l2.ldscore.gz"
    resources: slurm_extra = "--use-conda"
    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR_annotation_files/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{CHR}",
             snps = "../resources/ldsr/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for snRNAseq dwnSmpl: {wildcards.CELL_TYPE_GO_TERM}, CHR {wildcards.CHR}"
    log:     "../results/00LOG/LDSR/herring_GO_term_genes/ld_scores.snRNAseq.{CELL_TYPE_GO_TERM}.Chr{CHR}.log"
    shell:
        "python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule ldsr_stratified_baseline_v12_go_terms:
    input:   GWAS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/LDSR_annotation_files/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{CHR}.l2.ldscore.gz", CELL_TYPE_GO_TERM = config["CELL_TYPE_GO_TERM"], CHR = range(1,23))
    output:  "../results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{GWAS}_baseline.v1.2.results"
    resources: slurm_extra = "--use-conda"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR_annotation_files/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.",
             out_file = "../results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE_GO_TERM} and {wildcards.GWAS} GWAS"
    log:     "../results/00LOG/LDSR/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"
             
rule ldsr_stratified_summary_go_terms:
    input:   LDSR = expand("../results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/snRNAseq.{CELL_TYPE_GO_TERM}.{GWAS}_baseline.v1.2.results", CELL_TYPE_GO_TERM = config["CELL_TYPE_GO_TERM"], GWAS = config["GWAS"])
    output:  "../results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    resources: slurm_extra = "--use-conda"
    message: "Creating summary file for snRNAseq dwnSmpl: {wildcards.GWAS} GWAS"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/",
             cell_types_go_terms = "../resources/sheets/celltypes_goterms.tsv"
    log:     "../results/00LOG/LDSR/herring_GO_term_genes/snRNAseq.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """


             head -1 {params.dir}snRNAseq.L4_RORB_LRRK1-GO.0050877.lvl2.100UP_100DOWN.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types_go_terms}
              Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 ../results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/snRNAseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output} 2> {log}
             done

             """

