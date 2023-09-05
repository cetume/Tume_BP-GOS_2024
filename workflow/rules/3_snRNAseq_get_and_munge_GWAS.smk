localrules: get_and_munge_GWAS_for_MAGMA

## Download GWAS data
## SCZ data wave 3 v2 file is here: https://doi.org/10.6084/m9.figshare.14672178
## Height data Wood et al + UK BioBank is here: https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
## Munge GWAS data for MAGMA and SLDSR

rule get_and_munge_GWAS_for_MAGMA:
    output: SCZ = "../results/GWAS_for_MAGMA/SCZ_hg19_magma_ready.tsv",
            HEIGHT = "../results/GWAS_for_MAGMA/HEIGHT_hg19_magma_ready.tsv" 
    params: outdir = "../results/GWAS_for_MAGMA/"
    log:     "../results/logs/herring/snRNAseq.download.GWAS.log"
    shell:
             """
             #Download GWAS
             wget https://figshare.com/ndownloader/files/28169757 -P {params.outdir} #SCZ
             wget https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz -P {params.outdir} #HEIGHT
             
             #Unpack GWAS
             mv {params.outdir}28169757 {params.outdir}PGC3_SCZ_wave3_public.v2.tsv.gz
             gunzip {params.outdir}PGC3_SCZ_wave3_public.v2.tsv.gz
             gunzip {params.outdir}Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz
             
             #Remove SNP rs148878475
             sed -i '/rs148878475/d' {params.outdir}PGC3_SCZ_wave3_public.v2.tsv
             #Add N column (n = 161405)
             awk '{{s=(NR==1)?"N":"161405";$0=$0 OFS s}}1' {params.outdir}PGC3_SCZ_wave3_public.v2.tsv > {params.outdir}PGC3_SCZ_wave3_public.v2_2.tsv
             
             #Take correct columns and give new name
             awk '{{print $2"\t"$1"\t"$3"\t"$11"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$20}}' {params.outdir}PGC3_SCZ_wave3_public.v2_2.tsv|\
             sed 's/POS/BP/g' > {params.outdir}SCZ_hg19_magma_ready.tsv
             awk '{{print $3"\t"$1"\t"$2"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10}}' {params.outdir}Meta-analysis_Wood_et_al+UKBiobank_2018.txt |\
             sed 's/POS/BP/g' > {params.outdir}HEIGHT_hg19_magma_ready.tsv
             sed -i -e '1s/Tested_Allele/A1/' -e '1s/Other_Allele/A2/' {params.outdir}HEIGHT_hg19_magma_ready.tsv
             
             #Remove GWAS data
             rm {params.outdir}PGC3_SCZ_wave3_public.v2.tsv {params.outdir}PGC3_SCZ_wave3_public.v2_2.tsv {params.outdir}Meta-analysis_Wood_et_al+UKBiobank_2018.txt
             """

rule munge_GWAS_for_LDSR:
    input:  "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    output: "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz"
    params: merge = "../resources/ldsr/reference_files/w_hm3.snplist",
            out = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready"   
    resources: slurm_extra = "--use-conda"
    conda:   "../envs/ldsr.yml"
    log:     "../results/logs/herring/snRNAseq.munge.GWAS.{GWAS}.log"
    shell:
             """
             python ../resources/ldsr/munge_sumstats.py \
                    --sumstats {input} \
                    --merge-alleles {params.merge} \
                    --out {params.out} \
                    --a1-inc \
                    --chunksize 500000
             """
