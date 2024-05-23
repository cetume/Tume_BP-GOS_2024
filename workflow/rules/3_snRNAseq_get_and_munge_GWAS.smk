localrules: download_sumstats

rule download_sumstats:
    # Download GWAS sumstats files (see config)
    output:  "../results/GWAS/{GWAS}_hg19_raw.tsv"
    params:  lambda wildcards: config['GWAS'][wildcards.GWAS]
    message: "Download sumstats file"
    log:     "../results/00LOG/get_and_munge_GWAS/{GWAS}_download_sumstats.log"
    run:

             if wildcards.GWAS in ("MDD", "NEUROTICISM"):

                 shell("""

                 cp {params} temp; gunzip -c temp > {output}; rm temp 2> {log}

                 """)
             
             elif wildcards.GWAS in ("SCZ_EUR_ONLY", "BPD"):

                 shell("""

                 wget -O - {params} | gunzip -c | sed '/##/d' > {output} 

                 """)

             else:
                 
                 shell("""

                 wget -O - {params} | gunzip -c > {output} 

                 """)

rule standardise_sumstats:
    # Standardises sumstats: SNP, CHR. BP, PVAL, A1, A2 + additional GWAS dependant cols
    # python convert available here: https://github.com/precimed/python_convert/tree/master
    input:   "../results/GWAS/{GWAS}_hg19_raw.tsv"
    output:  "../results/GWAS/{GWAS}_hg19_standardised.tsv"
    message: "Formatting {input}"
    params: "../results/GWAS/{GWAS}_hg19_raw_temp.tsv"
    log:     "../results/00LOG/get_and_munge_GWAS/{GWAS}_standardise_sumstats.log"
    run:
             if wildcards.GWAS in ("HEIGHT", "BMI"):

                 shell("""
                 cat {input} | sed 's/Tested_Allele/A1/g' | sed 's/Other_Allele/A2/g' > {params};
    
                 python ../resources/python_convert/sumstats.py csv \
                 --sumstats {params} \
       	         --out {output} --force --auto --head 5 \
                 --log {log};

                 rm {params}
             
                  """)

             elif wildcards.GWAS in ("BPD", "SCZ_EUR_ONLY"):

                 shell("""
                 cat {input} | sed 's/ID/SNP/g' | sed 's/#CHROM/CHR/g' > {params};

                 python ../resources/python_convert/sumstats.py csv \
                 --sumstats {params} \
                 --out {output} --force --auto --head 5 \
                 --log {log};

                 rm {params}

                  """)             

             elif "SCZ" in wildcards.GWAS:

                 shell("""
                 sed -i '/rs148878475/d' {input};

                 python ../resources/python_convert/sumstats.py csv \
                 --sumstats {input} \
                 --out {output} --force --auto --head 5 \
                 --log {log};

                  """)

             elif "NEUROTICISM" in wildcards.GWAS:

                 shell("""
 
                 cat {input} | sed 's/POS/BP/g' | sed 's/RSID_UKB/SNP/g' | sed 's/REF/A1/g' | sed 's/ALT/A2/g' > {params};

                 python ../resources/python_convert/sumstats.py csv \
                 --sumstats {params} \
                 --out {output} --force --auto --head 5 \
                 --log {log};

                 rm {params}

                 """)

             else:

                 shell("""

                 python ../resources/python_convert/sumstats.py csv \
                 --sumstats {input} \
                 --out {output} --force --auto --head 5 \
                 --log {log}

                  """)

rule add_z_score:
    # Adds z-scores to GWAS sumstats lacking
    input:   "../results/GWAS/{GWAS}_hg19_standardised.tsv"
    output:  "../results/GWAS/{GWAS}_hg19_withZ.tsv"
    message: "Adding Z score to {input} if required"
    log:     "../results/00LOG/get_and_munge_GWAS/{GWAS}_addZscore.log"
    shell:
             """
             python ../resources/python_convert/sumstats.py zscore \
             --sumstats {input} \
             --out {output} --force \
             --log {log} \
             --a1-inc
             """

rule add_N:
    # N to GWAS sumstats lacking 
    input:   "../results/GWAS/{GWAS}_hg19_withZ.tsv"
    output:  "../results/GWAS/{GWAS}_hg19_withN.tsv"
    message: "Adding N to {input} if required"
    log:     "../results/00LOG/get_and_munge_GWAS/{GWAS}_addN.log"
    run:

             if "SCZ_EUR_ONLY" in wildcards.GWAS:

                 shell("""

                 awk '{{s=(NR==1)?"N":"130644";$0=$0 OFS s}}1' {input} > {output} 2> {log}

                 """)

             elif "SCZ" in wildcards.GWAS:
             
                 shell("""

                 awk '{{s=(NR==1)?"N":"161405";$0=$0 OFS s}}1' {input} > {output} 2> {log}
 
                 """)

             elif "ADHD" in wildcards.GWAS:

                 shell("""

                 awk '{{s=(NR==1)?"N":"225534";$0=$0 OFS s}}1' {input} > {output} 2> {log}
 
                 """)

             elif "ASD" in wildcards.GWAS:

                 shell("""

                 awk '{{s=(NR==1)?"N":"46350";$0=$0 OFS s}}1' {input} > {output} 2> {log}
 
                 """)
             
             elif "BPD" in wildcards.GWAS:

                 shell("""

                 awk '{{s=(NR==1)?"N":"413466";$0=$0 OFS s}}1' {input} > {output} 2> {log}

                 """)

             elif "NEUROTICISM" in wildcards.GWAS:

                 shell("""

                 awk '{{s=(NR==1)?"N":"313467";$0=$0 OFS s}}1' {input} > {output} 2> {log}

                 """)

             else:

                 shell("cp {input} {output}")

rule prep_for_magma:
    # Format sumstats for MAGMA input
    # Note that snakemake make objects (e.g. input) are available without curly braces here
    input:   "../results/GWAS/{GWAS}_hg19_withN.tsv"
    output:  "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    message: "Munging sumstats for MAGMA compatibility: {input}"
    log:     "../results/00LOG/GWAS/{GWAS}_prep_for_MAGMA.log"
    run:
             import pandas as pd
             import numpy as np
             df = pd.read_csv(input[0], sep = '\t')
             df = df.rename(columns={'PVAL': 'P'})
             first_cols = ['SNP','CHR','BP', 'P']
             last_cols = [col for col in df.columns if col not in first_cols]
             df = df[first_cols + last_cols]  
             df.to_csv(output[0], sep='\t', index = False, header = True)     

rule prep_for_ldsr:
    # Format sumstats for LDSR input
    input:   snps = "../resources/ldsr/reference_files/w_hm3.snplist",
             gwas = "../results/GWAS/{GWAS}_hg19_withN.tsv"
    output:  "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz"
    conda:   "../envs/ldsr.yml"
    message: "Munging sumstats for LDSR compatibility: {input.gwas}"
    params:  out = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready"
    log:     "../results/00LOG/GWAS/{GWAS}_hg19_prep_for_ldsr.log"
    shell:
        """

        python ../resources/ldsr/munge_sumstats.py --sumstats {input.gwas} \
        --merge-alleles {input.snps} \
        --out {params.out} \
        --a1-inc 2> {log}

        """
