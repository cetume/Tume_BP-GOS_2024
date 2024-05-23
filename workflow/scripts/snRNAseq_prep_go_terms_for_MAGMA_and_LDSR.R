#--------------------------------------------------------------------------------------
#
#    Prep GO terms enrichment test input files for MAGMA and SLDSR - L4-RORB-LRRK4
#
#--------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(R.utils)
library(readr)

## Set variables  ---------------------------------------------------------------------

SIG_CELLS <- 'L4_RORB_LRRK1'
GO_TERMS_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/resources/go_terms/'
DATA_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/gene_lists/herring/'
gene_coord <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds'
gene_coordinates <- readRDS(gene_coord)

## Load and filter GO Terms for L4-RORB-LRRK1 -----------------------------------------
# Top 20 significantly over-represented GO Terms

for (CELL_TYPE in SIG_CELLS) {
  
  CELL_TYPE_EDIT <- gsub("-", "_", CELL_TYPE)
  
  GO_TERMS <- read.csv(paste0(GO_TERMS_DIR, 'GO_terms_', CELL_TYPE, '.csv')) %>% 
    filter(Term == 'GO:0003008~system process' |
           Term == 'GO:0006813~potassium ion transport' |
           Term == 'GO:0007186~G-protein coupled receptor signaling pathway' |
           Term == 'GO:0007267~cell-cell signaling' |
           Term == 'GO:0007268~chemical synaptic transmission' |
           Term == 'GO:0050804~modulation of synaptic transmission' |
           Term == 'GO:0050808~synapse organization' |
           Term == 'GO:0071805~potassium ion transmembrane transport' |
           Term == 'GO:0050877~neurological system process' |
           Term == 'GO:0042391~regulation of membrane potential' |
           Term == 'GO:0006811~ion transport' |
           Term == 'GO:0098660~inorganic ion transmembrane transport' |
           Term == 'GO:0043269~regulation of ion transport' |
           Term == 'GO:0006836~neurotransmitter transport' |
           Term == 'GO:0019226~transmission of nerve impulse' |
           Term == 'GO:0055085~transmembrane transport' |
           Term == 'GO:0007270~neuron-neuron synaptic transmission' |
           Term == 'GO:0046903~secretion' |
           Term == 'GO:0007399~nervous system development' |
           Term == 'GO:0060078~regulation of postsynaptic membrane potential') %>%
    filter(FDR <= 0.05) %>%
    mutate(Term = gsub(' ', '_', Term)) %>%
    select(Term, Genes, Fold.Enrichment, FDR) 
  
  GO_TERMS_FILT <- GO_TERMS %>%
    select(Term, Genes) %>%
    group_by(Term) %>%
    dplyr::mutate(i1 = row_number()) %>%
    spread(Term, Genes) %>%
    select(-i1)

## Create enrichment files for MAGMA --------------------------------------------------
  
for (i in 1:ncol(GO_TERMS_FILT)) {
  
  TERM <- colnames(GO_TERMS_FILT)[i]
  GENES <- strsplit(as.character(GO_TERMS_FILT[1, i]), ", ")[[1]]
  
  if (file.exists(paste0(DATA_DIR, 'MAGMA/GO_term_genes_for_magma.txt'))) {
    
    cat('\n\nAPPEND\n\n')
    
    GENE_LIST <- paste0(CELL_TYPE_EDIT, '-', TERM, ' ', paste(GENES, collapse = ' '))
    cat('\n\n', GENE_LIST, '\n')
    cat(GENE_LIST, '\n', file = paste0(DATA_DIR, 'MAGMA/GO_term_genes_for_magma.txt'), 
        append = TRUE)
    
  } else {
    
    cat('\n\nCREATE FILE\n\n')
    
    GENE_LIST <- paste0(CELL_TYPE_EDIT, '-', TERM, ' ', paste(GENES, collapse = ' '))
    cat('\n', GENE_LIST, '\n')
    cat(GENE_LIST, '\n', file = paste0(DATA_DIR, 'MAGMA/GO_term_genes_for_magma.txt'))
    
  }
  
}

## Create enrichment files for SLDSR --------------------------------------------------

GO_TERMS_LDSR <- GO_TERMS %>%
  mutate(Term = paste0(CELL_TYPE_EDIT, '-', Term)) %>%
  separate_rows(Genes, convert = TRUE) %>%
  mutate(ensembl = Genes) %>%
  dplyr::select(Term, ensembl) 

GO_TERMS_LDSR <- GO_TERMS_LDSR %>%
  inner_join(gene_coordinates) %>% 
  mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
  dplyr::select(chr, start, end, ensembl, Term) %>%
  group_by(Term)
  
  GO_TERMS_LDSR$Term <- gsub(":", ".", GO_TERMS_LDSR$Term)
  GO_TERMS_LDSR$Term <- gsub("~.*", "", GO_TERMS_LDSR$Term)
  
  dir.create(paste0(DATA_DIR, 'LDSR_GO_term_genes/'),  recursive = TRUE, showWarnings = FALSE)
  GO_TERMS_LDSR %>% group_walk(~ write_tsv(.x[,1:4], paste0(DATA_DIR, 'LDSR_GO_term_genes/',
                                          .y$Term, '.lvl2.100UP_100DOWN.bed'), col_names = FALSE))

}
