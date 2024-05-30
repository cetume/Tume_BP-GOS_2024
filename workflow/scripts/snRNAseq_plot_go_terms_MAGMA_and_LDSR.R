#--------------------------------------------------------------------------------------
#
#    Plot MAGMA and SLDSR GO Term results
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(ggsignif)
library(cowplot)
library(reshape2)
library(ggplot2)
library(R.utils)

## Set variables  ---------------------------------------------------------------------

TRAIT <- 'SCZ_EUR_ONLY'

for(GWAS in TRAIT) {

magma <- paste0('~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma/snRNAseq_', GWAS, '.GO_term_genes.magma.35UP_10DOWN.gsa.out')
ldsr <- paste0('~/Desktop/Herring_snRNAseq_2023_pipeline/results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/snRNAseq_LDSR_', GWAS, '_baseline.v1.2_summary.tsv')
FIG_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/'
dir.create(paste0(FIG_DIR),  recursive = TRUE, showWarnings = FALSE)

## Load and prep data -----------------------------------------------------------------

  MAGMA_DF <- read.table(magma, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
  mutate(MAGMA = if_else(`BETA` > 0, MAGMA, 0)) %>%
  dplyr::select(FULL_NAME, MAGMA) %>%
  dplyr::rename(Category = FULL_NAME) %>%
  separate(Category, into=c('Category', 'GO_Term'), sep = '~', extra = "merge")
    
  MAGMA_DF$GO_Term<-capitalize(MAGMA_DF$GO_Term)
  
  MAGMA_DF <- MAGMA_DF %>%
    mutate(GO_Term = gsub('_', ' ', GO_Term))
  
  LDSR_FULL_DF <- read_tsv(ldsr) %>%
    mutate(SLDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0))
  
  LDSR_FULL_DF$Category <- sub("\\.", ":", LDSR_FULL_DF$Category)
  
  LDSR_FULL_DF <- LDSR_FULL_DF %>% filter(str_detect(Category, 'lvl2.100UP_100DOWN')) %>%
    separate(Category, into=c('Category', 'Window'), sep = '\\.', extra = "merge") 

  LDSR_DF <- LDSR_FULL_DF %>%
    dplyr::select(Category, SLDSR)

  PLOT_DF <- left_join(MAGMA_DF, LDSR_DF, by = 'Category') %>% 
    reshape2::melt() %>%
    separate(Category, into=c('Category', 'GO_code'), sep = '-', extra = "merge")
  
  PLOT_DF$Category <- gsub("_", "-", PLOT_DF$Category)

  levels <- c('Nervous system development',
              'System process',
              'Neurological system process',
              'G-protein coupled receptor signaling pathway',
              'Ion transport',
              'Potassium ion transport',
              'Transmembrane transport',
              'Potassium ion transmembrane transport',
              'Inorganic ion transmembrane transport',
              'Regulation of ion transport',
              'Synapse organization',
              'Neuron-neuron synaptic transmission',
              'Cell-cell signaling',
              'Chemical synaptic transmission',
              'Modulation of synaptic transmission',
              'Regulation of membrane potential',
              'Regulation of postsynaptic membrane potential',
              'Secretion',
              'Neurotransmitter transport', 
              'Transmission of nerve impulse')
  
## Create plots -----------------------------------------------------------------------
  
  BF_CORR <- 0.05/20

  MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(GO_Term, rev(levels)), width = 0.8,
                                                fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
    geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 16),
          axis.title.y = element_text(colour = "#000000", size = 16),
          axis.text.x = element_text(colour = "#000000", size = 13, vjust = 0.5),
          axis.text.y = element_text(colour = "#000000", size = 13),
          legend.text = element_text(colour = "#000000", size = 14),
          legend.position = "right",
          legend.title = element_blank(), 
          strip.text = element_text(size=14, face = 'bold')) +
    xlab(expression(-log[10](P))) +
    ylab('GO Terms') +
    xlim(0, 10) #+
      #facet_wrap('Category')
  
  assign(paste0(GWAS, '_magma_ldsr_GO_Term_lvl2_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv)
  
}

## Save plots --------------------------------------------------------------------------
  
  jpeg(file = paste0(FIG_DIR,'SCZ_EUR_ONLY_magma_ldsr_herring_GO_term_lvl2_plot.jpeg'), units = "in", width = 9, height = 8, res = 300)
  plot(SCZ_EUR_ONLY_magma_ldsr_GO_Term_lvl2_plot)
  dev.off()

