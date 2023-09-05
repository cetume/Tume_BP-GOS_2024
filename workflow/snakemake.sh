snakemake --profile ../config/profile/ $@ 2> smk-"`date +"%d-%m-%Y"`".log; mail -s "Snakemake has finished" tumece@cardiff.ac.uk < smk-"`date +"%d-%m-%Y"`".log
