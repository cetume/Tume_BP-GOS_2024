#!/bin/bash

mkdir -p ../results/GWAS_for_MAGMA
cd ../results/GWAS_for_MAGMA

#Download GWAS
wget https://figshare.com/ndownloader/files/28169757 #SCZ
wget https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz #HEIGHT


#Unpack GWAS
mv 28169757 PGC3_SCZ_wave3_public.v2.tsv.gz
gunzip PGC3_SCZ_wave3_public.v2.tsv.gz
gunzip Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz


#Remove SNP rs148878475
sed -i '/rs148878475/d' PGC3_SCZ_wave3_public.v2.tsv

#Add N column (n = 161405)
awk '{s=(NR==1)?"N":"161405";$0=$0 OFS s}1' PGC3_SCZ_wave3_public.v2.tsv > PGC3_SCZ_wave3_public.v2_2.tsv

#Take correct columns and give new name
awk '{print $2"\t"$1"\t"$3"\t"$11"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$20}' PGC3_SCZ_wave3_public.v2_2.tsv|\
sed 's/POS/BP/g' > SCZ_hg19_magma_ready.tsv

awk '{print $3"\t"$1"\t"$2"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10}' Meta-analysis_Wood_et_al+UKBiobank_2018.txt |\
sed 's/POS/BP/g' > HEIGHT_hg19_magma_ready.tsv

sed -i -e '1s/Tested_Allele/A1/ -e '1s/Other_Allele/A2/' HEIGHT_hg19_magma_ready.tsv

#Remove GWAS data
rm PGC3_SCZ_wave3_public.v2.tsv PGC3_SCZ_wave3_public.v2_2.tsv Meta-analysis_Wood_et_al+UKBiobank_2018.txt
