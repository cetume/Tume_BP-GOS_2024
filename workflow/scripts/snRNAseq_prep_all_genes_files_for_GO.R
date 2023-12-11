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
LEVEL <- 2

## Load and prepare data --------------------------------------------------------------
load(paste0(CTD_DIR, 'ctd_herring.rda')) 

CELL_TYPES <- colnames(ctd[[LEVEL]]$specificity_quantiles) 

# Generate tables for chosen cell populations
cell_types <- c("L4_RORB_dev-2", "L4_RORB_LRRK1") 

for (type in cell_types) { 
all_genes <- as_tibble(as.matrix(ctd[[LEVEL]]$specificity_quantiles), rownames = 'hgnc') %>%
  inner_join(gene_coordinates) %>%
  pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
  filter(quantile != 0) %>%
  filter(cell_type == type) %>%
  dplyr::select(ensembl) 

## Save txt files ---------------------------------------------------------------------
write.table(all_genes, paste0(OUT_DIR, type, '_all_genes.txt'), quote = 
FALSE, col.names = FALSE, row.names = FALSE)
}
