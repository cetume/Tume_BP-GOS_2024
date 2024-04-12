#--------------------------------------------------------------------------------------
#
#    Prepare all genes files for the background in GO Term analysis
#
#--------------------------------------------------------------------------------------

##Load packages  ----------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(EWCE)

##  Set variables  --------------------------------------------------------------------
CTD_DIR <- 'Herring_snRNAseq_2023_pipeline/results/ctd_object/' 
OUT_DIR <- 'Herring_snRNAseq_2023_pipeline/resources/go_terms/'
gene_coord <- 'Herring_snRNAseq_2023_pipeline/results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds'
LEVEL <- 2

## Load and prepare data --------------------------------------------------------------
load(paste0(CTD_DIR, 'ctd_herring.rda')) 

CELL_TYPES <- colnames(ctd[[LEVEL]]$specificity_quantiles) 

gene_coordinates <- readRDS(gene_coord)

# Generate tables for chosen cell populations
SIG_CELLS <- c("L4_RORB_dev-2", "L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_dev-1", "L4_RORB_dev-fetal", "L4_RORB_MME") 

for (TYPE in SIG_CELLS) { 
all_genes <- as_tibble(as.matrix(ctd[[LEVEL]]$specificity_quantiles), rownames = 'hgnc') %>%
  inner_join(gene_coordinates) %>%
  pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
  filter(quantile != 0) %>%
  filter(cell_type == TYPE) %>%
  dplyr::select(ensembl) 

## Save txt files ---------------------------------------------------------------------
write.table(all_genes, paste0(OUT_DIR, 'All_genes_expressed_', TYPE, '.txt'), quote = 
FALSE, col.names = FALSE, row.names = FALSE)
}
