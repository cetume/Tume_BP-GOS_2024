#--------------------------------------------------------------------------------------
#
#    Prep GO terms enrichment test input files for MAGMA and SLDSR
#
#--------------------------------------------------------------------------------------

#Script also plots the fold enrichment and FDR for the top 10 GO terms for each implicated cell population

## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(R.utils)

## Set variables  ---------------------------------------------------------------------

SIG_CELLS <- c('L4_RORB_LRRK1', 'L4_RORB_dev-2', 'L4_RORB_MET')
GO_TERMS_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/resources/go_terms/'
DATA_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/gene_lists/herring/'
FIG_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/'
gene_coord <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/R_objects/Ensembl_hg19_gene_coords_noMHC.rds'
gene_coordinates <- readRDS(gene_coord)

dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

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


## Plot GO term enrichment figure -----------------------------------------------------

# Prepare data
GO_TERMS_FOR_PLOT <- GO_TERMS %>% 
    separate(Term, into=c('GO_Code', 'GO_Term'), sep = '~', extra = "merge")
  
GO_TERMS_FOR_PLOT$GO_Term <- capitalize(GO_TERMS_FOR_PLOT$GO_Term)
  
GO_TERMS_FOR_PLOT <- GO_TERMS_FOR_PLOT %>%
    mutate(GO_Term = gsub('_', ' ', GO_Term)) %>%
    dplyr::select(GO_Term, Fold.Enrichment, FDR)
  
colnames(GO_TERMS_FOR_PLOT) <- c('Term', paste0(CELL_TYPE_EDIT, '-FE'), paste0(CELL_TYPE_EDIT, '-FDR'))
  
assign(paste0(CELL_TYPE_EDIT, '_GO_Terms'), GO_TERMS_FOR_PLOT, envir = .GlobalEnv)
  
}
  
df_list <- list(L4_RORB_dev_2_GO_Terms, L4_RORB_LRRK1_GO_Terms, L4_RORB_MET_GO_Terms)

GO_DATA <- df_list %>% reduce(full_join)

GO_DATA <- GO_DATA %>%
  pivot_longer(-Term) %>%
  separate(name, into = c("cell_type", "Score"), '-') %>%
  pivot_wider(names_from = Score, values_from = value)

colnames(GO_DATA) <- c("Term", "cell_type", "Fold Enrichment", "FDR")

GO_DATA$FDR <- as.numeric(as.character(GO_DATA$FDR))
GO_DATA$cell_type <- gsub("_", "-", GO_DATA$cell_type)

# Plot data 
levels <- c('Nervous system development', 
            'Neuron development', 
            'Cell projection organization',
            'System process', 
            'Neurological system process',
            'G-protein coupled receptor signaling pathway',
            'Metal ion transport', 
            'Potassium ion transport', 
            'Regulation of ion transport',
            'Regulation of ion transmembrane transport',
            #'Regulation of transmembrane transport',
            'Synapse organization',
            'Synapse assembly',
            #'Synaptic signaling',
            'Trans-synaptic signaling',
            'Cell-cell signaling',
            'Regulation of synapse structure or activity',
            'Modulation of synaptic transmission',
            'Regulation of synaptic plasticity',
            'Regulation of membrane potential', 
            'Detection of stimulus involved in sensory perception')

GO_PLOT <- ggplot(data = GO_DATA, aes(y = factor(Term, level = rev(levels)), x = cell_type,
                                      color = -log10(FDR), size = `Fold Enrichment`)) +
  geom_point() +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        #      panel.grid.major = element_blank(),
        #    panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 12, vjust = -0.5),
        axis.title.y = element_text(colour = "#000000", size = 12),
        axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.8, hjust = 0, angle = -45),
        axis.text.y  = element_text(colour = "#000000", size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(colour = "#000000")) +
  ylab("GO Terms") +
  xlab("Cell type") +
  scale_radius(limits = c(1, 4), range = c(1,10))

tiff(paste0(FIG_DIR, "Cell_pops_top10_GoTerm_enrichment.tiff"), height = 18, width = 19, units='cm',
     compression = "lzw", res = 300)
print(GO_PLOT)
dev.off()

