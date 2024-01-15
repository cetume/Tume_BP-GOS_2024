#--------------------------------------------------------------------------------------
#
#    Prepare enrichment file for MAGMA interaction analyses
#
#--------------------------------------------------------------------------------------

#Gene sets include the genes for 4 GO terms of interest (GO:0007399~nervous system development, GO:0042391~regulation of membrane potential, GO:0043269~regulation of ion transport, GO:0099537~trans-synaptic signaling) and genes in the top expression specificity decile for L4_RORB_LRRK1


##  Load Packages  --------------------------------------------------------------------

library(tidyverse)

## Set variables  ---------------------------------------------------------------------

DATA_DIR <- 'Herring_snRNAseq_2023_pipeline/resources/go_terms/'
CTD_DIR <- 'Herring_snRNAseq_2023_pipeline/results/ctd_objects/'
OUT_DIR <- 'Herring_snRNAseq_2023_pipeline/results/gene_lists/herring/MAGMA/'
gene_coord <- 'Herring_snRNAseq_2023_pipeline/results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds'
gene_coordinates <- readRDS(gene_coord)
LEVEL <- 2

## Load and prep genes of GO terms ----------------------------------------------------

for (GO in c("GO:0099537", "GO:0007399", "GO:0043269", "GO:0042391")) {

GO_EDIT <- gsub(":", "", GO)

  GO_GENES <- read.csv(paste0(DATA_DIR, GO_EDIT, '_genes.csv')) %>%
  dplyr::rename(ensembl = "ID") %>%
  mutate(ensembl = ifelse(ensembl == "ENSG00000006837, ENSG00000273345", "ENSG00000006837",
         ifelse(ensembl == "ENSG00000251349, ENSG00000241697", "ENSG00000241697", ensembl))) %>%
  inner_join(gene_coordinates) %>%
  dplyr::select(ensembl) %>%
  filter(ensembl != "") %>%
  mutate(term = paste0(GO))
  
  assign(paste0(GO_EDIT, '_Genes'), GO_GENES, envir = .GlobalEnv)  
  
}

## Load and prep genes of L4_RORB_LRRK1 -----------------------------------------------

load(paste0(CTD_DIR, 'ctd_herring.rda'))

CELL_TYPES <- colnames(ctd[[LEVEL]]$specificity_quantiles)

LRRK1_Genes <- as_tibble(as.matrix(ctd[[LEVEL]]$specificity_quantiles), rownames = 'hgnc') %>%
  inner_join(gene_coordinates) %>%
  pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
  filter(quantile == 10) %>%
  filter(cell_type == "L4_RORB_LRRK1") %>%
  dplyr::select(cell_type, ensembl)
  
colnames(LRRK1_Genes) <- c('term', 'ensembl')

## Merge data frames and generate MAGMA enrichment file -------------------------------

GO_GENE_DATA <- merge(GO0099537_Genes, GO0007399_Genes_test, all = TRUE) %>%
  merge(GO0043269_Genes, all = TRUE) %>%
  merge(GO0042391_Genes, all = TRUE) %>%
  merge(LRRK1_Genes, all = TRUE) %>%
  with(., split(ensembl, term))

for(i in names(GO_GENE_DATA)) {
  
  cat(i, " ", paste(GO_GENE_DATA[[i]], collapse = " "), "\n",
      file = paste0(OUT_DIR, 'GO_term_genes_for_cond_magma.txt'), sep = '', append = TRUE)
  
}

