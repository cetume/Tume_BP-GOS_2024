#!/usr/bin/env bash
#SBATCH --job-name=install_packages
#SBATCH --partition highmem
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node 1
#SBATCH --account=scw2043
#SBATCH --output /scratch/c.c1837163/logs/%J.o
#SBATCH --error /scratch/c.c1837163/logs/%J.e
#SBATCH --mail-type=END
#SBATCH --mail-user=tumece@cardiff.ac.uk


module load R/3.5.1

Rscript --vanilla /scratch/c.c1837163/Herring_snRNAseq_2023/workflow/scripts/snRNAseq_generate_ctd.R
