#--------------------------------------------------------------------------------------
#
#    Plot MAGMA GO Term results
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------

library(dplyr)
library(ggsignif)
library(cowplot)
library(reshape2)
library(ggplot2)
library(R.utils)

## Set variables  ---------------------------------------------------------------------

magma <- 'Herring_snRNAseq_2023/results/magma/snRNAseq_SCZ.GO_term_genes.magma.35UP_10DOWN.gsa.out')
FIG_DIR <- 'Herring_snRNAseq_2023/results/figures/'

## Load and prep data -----------------------------------------------------------------

  MAGMA_DF <- read.table(magma, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
  dplyr::select(FULL_NAME, MAGMA) %>%
  dplyr::rename(Category = FULL_NAME) %>%
  separate(Category, into=c('Category', 'GO_Term'), sep = '-', extra = "merge") %>%
  separate(GO_Term, into=c('GO_code', 'GO_Term'), sep = '~', extra = "merge")
    
    MAGMA_DF$GO_Term<-capitalize(MAGMA_DF$GO_Term)
  
  MAGMA_DF <- MAGMA_DF %>%
  mutate(GO_Term = gsub('_', ' ', GO_Term))
  
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
              'Regulation of transmembrane transport',
              'Synapse organization',
              'Synapse assembly',
              'Synaptic signaling',
              'Trans-synaptic signaling',
              'Cell-cell signaling',
              'Regulation of synapse structure or activity',
              'Modulation of synaptic transmission',
              'Regulation of synaptic plasticity',
              'Regulation of membrane potential', 
              'Detection of stimulus involved in sensory perception')
  
## Creating plot ----------------------------------------------------------------------
  
  MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(x = MAGMA, y = factor(GO_Term, rev(levels)),
                                            fill = '#F8766D')) +
  geom_bar(stat = "identity", color = 'black', position = "dodge") +
  #geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
 # ggtitle(paste0('GO Terms')) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 15),
        axis.title.y = element_text(colour = "#000000", size = 15),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.position = "none",
        strip.text = element_text(size=14, face = 'bold')) +
  xlab(expression(-log[10](P))) +
  ylab('GO Terms') +
  xlim(0, 13) +
    facet_wrap('Category')
  
## Save plot --------------------------------------------------------------------------
  
  jpeg(file = paste0(FIG_DIR,'SCZ_magma_herring_GO_term_lvl2_plot.jpeg'), units = "in", width = 12, height = 8, res = 300)
  plot(MAGMA_PLOT)
  dev.off()
