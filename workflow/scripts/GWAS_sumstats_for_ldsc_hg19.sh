#!/bin/bash

python ../resources/ldsr/munge_sumstats.py \
        --sumstats ../results/GWAS_for_MAGMA/HEIGHT_hg19_magma_ready.tsv \
        --merge-alleles ../resources/ldsr/reference_files/w_hm3.snplist \
        --out ../results/GWAS_for_ldsc/HEIGHT_hg19_ldsc_ready \
        --a1-inc \
        --chunksize 500000
