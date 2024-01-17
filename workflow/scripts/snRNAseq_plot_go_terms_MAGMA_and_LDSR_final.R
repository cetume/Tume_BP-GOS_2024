#--------------------------------------------------------------------------------------
#
#    Plot MAGMA GO Term results
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

magma <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma/snRNAseq_SCZ.GO_term_genes.magma.35UP_10DOWN.gsa.out'
ldsr <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/LDSR_part_herit/baseline_v1.2/herring_GO_term_genes/snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv'
FIG_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/'
dir.create(paste0(FIG_DIR),  recursive = TRUE, showWarnings = FALSE)

## Load and prep data -----------------------------------------------------------------

  MAGMA_DF <- read.table(magma, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
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
    separate(Category, into=c('Category', 'GO_code'), sep = '-', extra = "merge") %>%
    filter(GO_Term != 'Regulation of transmembrane transport' & GO_Term != 'Synaptic signaling')
  
  PLOT_DF$Category <- gsub("_", "-", PLOT_DF$Category)

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
  
## Creating plots ---------------------------------------------------------------------
  
  BF_CORR <- 0.05/50

  MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(GO_Term, rev(levels)),
                                                fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
    #geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    #ggtitle(GWAS) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 15),
          axis.title.y = element_text(colour = "#000000", size = 15),
          axis.text.x = element_text(colour = "#000000", size = 12, vjust = 0.5),
          axis.text.y = element_text(colour = "#000000", size = 12),
          legend.position = "right",
          legend.title = element_blank(), 
          strip.text = element_text(size=14, face = 'bold')) +
    xlab(expression(-log[10](P))) +
    ylab('GO Terms') +
    xlim(0, 13) +
      facet_wrap('Category')

  PLOT_mean <- PLOT_DF %>% pivot_wider(names_from = variable, values_from = value)
  PLOT_mean$mean <- rowMeans(PLOT_mean[,c('MAGMA', 'SLDSR')])
  PLOT_mean <- PLOT_mean %>% mutate(COLOUR = ifelse(MAGMA > -log10(BF_CORR) & SLDSR > -log10(BF_CORR), "Both",
                                             ifelse(MAGMA > -log10(BF_CORR), "MAGMA",
                                             ifelse(SLDSR > -log10(BF_CORR), "SLDSR", "None"))))
  
  colour_table <- tibble(
    COLOUR = c("Both", "MAGMA", "SLDSR", "None"),
    Code = c("#00BA38", "yellow", "#00B0F6", "lightgrey")
  )

  PLOT_mean$COLOUR <- factor(PLOT_mean$COLOUR, levels = colour_table$COLOUR)

  MAGMA_LDSR_MEAN_PLOT <- ggplot(data = PLOT_mean, aes(x = mean, y = factor(GO_Term, rev(levels)), fill = COLOUR)) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    scale_fill_manual(values = colour_table$Code, drop = FALSE, name = "Significant at 
Bonferroni threshold") +
    geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
    #geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    #ggtitle(GWAS) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 15),
          axis.title.y = element_text(colour = "#000000", size = 15),
          axis.text.x = element_text(colour = "#000000", size = 12, vjust = 0.5),
          axis.text.y = element_text(colour = "#000000", size = 14),
          title = element_text(colour = "#000000", size = 16),
          legend.text = element_text(colour = "#000000", size = 14),
          legend.position = "right",
          strip.text = element_text(size=14, face = 'bold')) +
    xlab(expression(Mean -log[10](P))) +
    ylab('GO Terms') +
    xlim(0, 8.5) +
    facet_wrap('Category')

## Save plots --------------------------------------------------------------------------

  #jpeg(file = paste0(FIG_DIR,'SCZ_magma_ldsr_herring_GO_term_lvl2_plot.jpeg'), units = "in", width = 12, height = 8, res = 300)
  #plot(MAGMA_LDSR_PLOT)
  #dev.off()
  
  jpeg(file = paste0(FIG_DIR,'SCZ_magma_ldsr_mean_herring_GO_term_lvl2_plot.jpeg'), units = "in", width = 15, height = 8, res = 300)
  plot(MAGMA_LDSR_MEAN_PLOT)
  dev.off()

