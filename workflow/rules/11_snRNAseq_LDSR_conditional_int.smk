rule ldsr_stratified_baseline_v12_internal_cond:
    input:   GWAS = "../results/GWAS_for_ldsc/SCZ_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/LDSR_annotation_files/herring/snRNAseq.{COND_CELL_A}.{CHR}.l2.ldscore.gz", COND_CELL_A = config["COND_CELL_A"], CHR = range(1,23)),
             COND = expand("../results/LDSR_annotation_files/herring/snRNAseq.{COND_CELL_B}.{CHR}.l2.ldscore.gz", COND_CELL_B = config["COND_CELL_B"], CHR = range(1,23))
    output:  "../results/LDSR_part_herit/baseline_v1.2/herring_cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.SCZ_baseline.v1.2.results"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR_annotation_files/herring/snRNAseq.{COND_CELL_A}.",
             cond_anns = "../results/LDSR_annotation_files/herring/snRNAseq.{COND_CELL_B}.",
             out_file = "../results/LDSR_part_herit/baseline_v1.2/herring_cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.SCZ_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.COND_CELL_A} vs. {wildcards.COND_CELL_B}, 100UP_100DOWN and SCZ"
    log:     "../results/00LOG/LDSR/herring_cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.SCZ.baseline.v1.2_partHerit.log"
    shell:
             """
             if [[ {wildcards.COND_CELL_A} == {wildcards.COND_CELL_B} ]] 
             then 
               printf 'L2_2\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n' > {output} 2> {log} 
             else  
               python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} --ref-ld-chr {params.baseline},{params.cond_anns},{params.LD_anns} --overlap-annot --frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}
             fi
             """ 

rule ldsr_stratified_summary_internal_cond:
    input:   expand("../results/LDSR_part_herit/baseline_v1.2/herring_cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.SCZ_baseline.v1.2.results", COND_CELL_A = config['COND_CELL_A'], COND_CELL_B = config['COND_CELL_B'])
    output:  "../results/LDSR_part_herit/baseline_v1.2/herring_cond_int/snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv"
    message: "Creating internal conditional summary file for SCZ"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/herring_cond_int/",
             cell_types = "../resources/sheets/celltypes_conditional_int.tsv"
    log:     "../results/00LOG/LDSR/herring_cond_int/snRNAseq.SCZ_baseline.v1.2_partHerit.summary.log"
    shell:
             """


             head -1 {params.dir}snRNAseq.L4_RORB_dev-2.lvl2.100UP_100DOWN_vs_L4_RORB_LRRK1.lvl2.100UP_100DOWN.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_2 ../results/LDSR_part_herit/baseline_v1.2/herring_cond_int/snRNAseq."$Line".SCZ_baseline.v1.2.results | sed "s/L2_2/$Line/g" >> {output} 2> {log}
             done

             """
