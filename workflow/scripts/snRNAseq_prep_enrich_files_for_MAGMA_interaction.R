#--------------------------------------------------------------------------------------
#
#    Prepare enrichment file for MAGMA interaction analyses
#
#--------------------------------------------------------------------------------------

#Gene sets include all genes from the 10 enriched GO terms and genes in the top expression specificity decile for L4_RORB_LRRK1
#GO terms includes are: GO:0003008~system process, GO:0050877~neurological system process, GO:0007267~cell-cell signaling, GO:0099537~trans-synaptic signaling, GO:0030001~metal ion transport, GO:0007399~nervous system development, GO:0048666~neuron development, GO:0034765~regulation of ion transmembrane transport, GO:0043269~regulation of ion transport, GO:0042391~regulation of membrane potential

##  Load Packages  --------------------------------------------------------------------

library(tidyverse)

## Set variables  ---------------------------------------------------------------------

DATA_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/resources/go_terms/'
CTD_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/ctd_objects/'
OUT_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/gene_lists/herring/MAGMA/'
gene_coord <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds'
gene_coordinates <- readRDS(gene_coord)
LEVEL <- 2

## Load and prep genes of GO terms ----------------------------------------------------

for (GO in c("GO:0003008", "GO:0050877", "GO:0007267", "GO:0099537", "GO:0030001", "GO:0007399", "GO:0048666", "GO:0034765", "GO:0043269", "GO:0042391")) {

GO_EDIT <- gsub(":", "", GO)

  GO_GENES <- read.csv(paste0(DATA_DIR, GO_EDIT, '_genes.csv')) %>%
  dplyr::rename(ensembl = "ID") %>%
  inner_join(gene_coordinates) %>%
  dplyr::select(ensembl) %>%
  filter(ensembl != "") %>%
  mutate(Category = paste0(GO))
  
  assign(paste0(GO_EDIT), GO_GENES, envir = .GlobalEnv)  
  
}

## Load and prep genes of L4_RORB_LRRK1 -----------------------------------------------

load(paste0(CTD_DIR, 'ctd_herring.rda'))

CELL_TYPES <- colnames(ctd[[LEVEL]]$specificity_quantiles)

LRRK1 <- as_tibble(as.matrix(ctd[[LEVEL]]$specificity_quantiles), rownames = 'hgnc') %>%
  inner_join(gene_coordinates) %>%
  pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
  filter(quantile == 10) %>%
  filter(cell_type == "L4_RORB_LRRK1") %>%
  dplyr::select(cell_type, ensembl)
  
colnames(LRRK1) <- c('Category', 'ensembl')

## Merge data frames and generate MAGMA enrichment file -------------------------------

df_list <- list(GO0003008, GO0050877, GO0007267, GO0099537, GO0030001, GO0007399, GO0048666, GO0034765, GO0043269, GO0042391, LRRK1)

GO_GENE_DATA <- df_list %>% reduce(full_join) %>%
  with(., split(ensembl, Category))

for(i in names(GO_GENE_DATA)) {
  
  cat(i, " ", paste(GO_GENE_DATA[[i]], collapse = " "), "\n",
      file = paste0(OUT_DIR, 'GO_term_genes_for_cond_magma_all.txt'), sep = '', append = TRUE)
  
}

