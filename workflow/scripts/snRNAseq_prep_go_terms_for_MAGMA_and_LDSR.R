#--------------------------------------------------------------------------------------
#
#    Prep GO terms enrichment test input files for MAGMA and SLDSR
#
#--------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------

library(dplyr)
library(tidyr)

## Set variables  ---------------------------------------------------------------------

SIG_CELLS <- c('L4_RORB_LRRK1', 'L4_RORB_dev-2', 'L4_RORB_MET')
GO_TERMS_DIR <- 'Herring_snRNAseq_2023_pipeline/resources/go_terms/'
DATA_DIR <- 'Herring_snRNAseq_2023_pipeline/results/gene_lists/herring/'
gene_coord <- 'Herring_snRNAseq_2023_pipeline/results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds'
gene_coordinates <- readRDS(gene_coord)

## Create enrichment files for MAGMA --------------------------------------------------

for (CELL_TYPE in SIG_CELLS) {
  
CELL_TYPE_EDIT <- gsub("-", "_", CELL_TYPE)

GO_TERMS <- read.csv(paste0(GO_TERMS_DIR, 'GO_terms_', CELL_TYPE, '.csv')) %>% filter(Term == 'GO:0099537~trans-synaptic signaling' |
                            Term == 'GO:0050804~modulation of synaptic transmission' |
                            Term == 'GO:0007267~cell-cell signaling' |
                            Term == 'GO:0043269~regulation of ion transport' |
                            Term == 'GO:0050808~synapse organization' |
                            Term == 'GO:0048167~regulation of synaptic plasticity' |
                            Term == 'GO:0003008~system process' |
                            Term == 'GO:0006813~potassium ion transport' |
                            Term == 'GO:0050877~neurological system process' |
                            Term == 'GO:0007186~G-protein coupled receptor signaling pathway' |
                            Term == 'GO:0042391~regulation of membrane potential' |
                            Term == 'GO:0034765~regulation of ion transmembrane transport'|
                            Term == 'GO:0030001~metal ion transport' |
                            Term == 'GO:0007399~nervous system development' |
                            Term == 'GO:0048666~neuron development' |
                            Term == 'GO:0030030~cell projection organization' |
                            Term == 'GO:0050803~regulation of synapse structure or activity' |
                            Term == 'GO:0007416~synapse assembly' |
                            #Term == 'GO:0099536~synaptic signaling' |
                            #Term == 'GO:0034762~regulation of transmembrane transport' |
                            Term == 'GO:0050906~detection of stimulus involved in sensory perception') %>%
  filter(FDR <= 0.05) %>%
  mutate(Term = gsub(' ', '_', Term)) %>%
  select(Term, Genes, Fold.Enrichment, FDR)

GO_TERMS_FILT <- GO_TERMS %>%
  select(Term, Genes) %>%
  group_by(Term) %>%
  dplyr::mutate(i1 = row_number()) %>%
  spread(Term, Genes) %>%
  select(-i1)

for (i in 1:ncol(GO_TERMS_FILT)) {
  
  TERM <- colnames(GO_TERMS_FILT)[i]
  GENES <- strsplit(as.character(GO_TERMS_FILT[1, i]), ", ")[[1]]
  
  if (file.exists(paste0(DATA_DIR, 'GO_term_genes_for_magma.txt'))) {
    
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

GO_TERMS_2 <- GO_TERMS %>%
  mutate(Term = paste0(CELL_TYPE_EDIT, '-', Term)) %>%
  separate_rows(Genes, convert = TRUE) %>%
  dplyr::select(Term, Genes) colnames(GO_TERMS_2) <- c("cell_type", "ensembl")

GO_TERMS_2 <- GO_TERMS_2 %>%
  inner_join(gene_coordinates) %>%
  mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
  dplyr::select(chr, start, end, ensembl, cell_type) %>%
  group_by(cell_type)
  
  GO_TERMS_2$cell_type <- gsub(":", ".", GO_TERMS_2$cell_type)
  GO_TERMS_2$cell_type <- gsub("~.*", "", GO_TERMS_2$cell_type)
  GO_TERMS_2 %>% group_walk(~ write_tsv(.x[,1:4], paste0(DATA_DIR, 'LDSR_GO_term_genes/',
                                          .y$cell_type, '.lvl2.100UP_100DOWN.bed'), col_names = FALSE))
}