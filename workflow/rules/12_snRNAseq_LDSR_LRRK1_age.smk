localrules: ldsr_stratified_summary_LRRK1_age

rule ldsr_make_annot_LRRK1_age:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   gene_set = "../results/gene_lists/herring/LDSR/{LRRK1_AGE}.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/LDSR_annotation_files/herring/snRNAseq.{LRRK1_AGE}.{CHR}.annot.gz"
    resources: slurm_extra = "--use-conda"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snRNAseq: {wildcards.LRRK1_AGE}, Chr {wildcards.CHR}"
    log:     "../results/00LOG/LDSR/herring/make_annot.snRNAseq.{LRRK1_AGE}.Chr{CHR}.log"
    shell:
             """

             python ../resources/ldsr/make_annot.py \
             --bed-file {input.gene_set} \
             --windowsize 0 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}

             """

rule ldsr_ld_scores_LRRK1_age:
    input:   annot = "../results/LDSR_annotation_files/herring/snRNAseq.{LRRK1_AGE}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsr/reference_files/hapmap3_snps"
    output:  "../results/LDSR_annotation_files/herring/snRNAseq.{LRRK1_AGE}.{CHR}.l2.ldscore.gz"
    resources: slurm_extra = "--use-conda"
    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR_annotation_files/herring/snRNAseq.{LRRK1_AGE}.{CHR}",
             snps = "../resources/ldsr/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for snRNAseq dwnSmpl: {wildcards.LRRK1_AGE}, CHR {wildcards.CHR}"
    log:     "../results/00LOG/LDSR/herring/ld_scores.snRNAseq.{LRRK1_AGE}.Chr{CHR}.log"
    shell:
        "python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"


rule ldsr_stratified_baseline_v12_LRRK1_age:
    input:   GWAS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/LDSR_annotation_files/herring/snRNAseq.{LRRK1_AGE}.{CHR}.l2.ldscore.gz", LRRK1_AGE = config["LRRK1_AGE"], CHR = range(1,23))
    output:  "../results/LDSR_part_herit/baseline_v1.2/herring/snRNAseq.{LRRK1_AGE}.{GWAS}_baseline.v1.2.results"
    resources: slurm_extra = "--use-conda"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR_annotation_files/herring/snRNAseq.{LRRK1_AGE}.",
             out_file = "../results/LDSR_part_herit/baseline_v1.2/herring/snRNAseq.{LRRK1_AGE}.{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.LRRK1_AGE} and {wildcards.GWAS} GWAS"
    log:     "../results/00LOG/LDSR/herring/snRNAseq.{LRRK1_AGE}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule ldsr_stratified_summary_LRRK1_age:
    input:   LDSR = expand("../results/LDSR_part_herit/baseline_v1.2/herring/snRNAseq.{LRRK1_AGE}.{GWAS}_baseline.v1.2.results", LRRK1_AGE = config["LRRK1_AGE"], GWAS = config["GWAS"])
    output:  "../results/LDSR_part_herit/baseline_v1.2/herring/snRNAseq_LDSR_LRRK1_age_{GWAS}_baseline.v1.2_summary.tsv"
    resources: slurm_extra = "--use-conda"
    message: "Creating summary file for snRNAseq: {wildcards.GWAS} GWAS"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/herring/",
             ages = "../resources/sheets/LRRK1_age.tsv"
    log:     "../results/00LOG/LDSR/herring/snRNAseq.LRRK1_age_{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """


             head -1 {params.dir}snRNAseq.Infancy.100UP_100DOWN.SCZ_EUR_ONLY_baseline.v1.2.results > {output}
             File={params.ages}
              Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 ../results/LDSR_part_herit/baseline_v1.2/herring/snRNAseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output} 2> {log}
             done

             """

