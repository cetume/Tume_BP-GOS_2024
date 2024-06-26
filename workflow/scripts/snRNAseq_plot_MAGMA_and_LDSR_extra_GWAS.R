#--------------------------------------------------------------------------------------
#
#    Plot MAGMA and SLDSR results - Additional GWAS traits
#
#--------------------------------------------------------------------------------------

## Load Packages ----------------------------------------------------------------------
library(tidyverse)
library(ggsignif)
library(cowplot)
library(reshape2)
library(Seurat)
library(ggforce)

## Set variables ----------------------------------------------------------------------
cat('\nPlotting MAGMA and SLDSR results ... \n')


TRAIT <- c('ADHD', 'ASD', 'BPD', 'MDD', 'NEUROTICISM', 'HEIGHT')
PROJECTS <- 'herring'
level <- '2'
FIG_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/'
dir.create(paste0(FIG_DIR),  recursive = TRUE, showWarnings = FALSE)

BF_CORR <- 0.05/84

for(PROJECT in PROJECTS) {
  for(GWAS in TRAIT) {
    
## Load data --------------------------------------------------------------------------
    
magma <- paste0('~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma/snRNAseq_', GWAS, '.', PROJECT, '.lvl', level, '.magma.35UP_10DOWN.gsa.out')
ldsr <- paste0('~/Desktop/Herring_snRNAseq_2023_pipeline/results/LDSR_part_herit/baseline_v1.2/', PROJECT, '/snRNAseq_LDSR_', GWAS, '_baseline.v1.2_summary.tsv')
    
## Preparing data ---------------------------------------------------------------------
cat('\nPreparing data ... \n')
    
MAGMA_DF <- read.table(magma, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
  mutate(MAGMA = if_else(`BETA` > 0, MAGMA, 0)) %>%
  dplyr::select(VARIABLE, MAGMA) %>%
  dplyr::rename(Category = VARIABLE)
    
LDSR_FULL_DF <- read_tsv(ldsr) %>%
  mutate(SLDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
  filter(str_detect(Category, paste0('lvl', level, '.100UP_100DOWN'))) %>%
  separate(Category, into=c('Category', 'Window'), sep = '\\.', extra = "merge")
    
LDSR_DF <- LDSR_FULL_DF %>%
  dplyr::select(Category, SLDSR)
    
## Creating plots ---------------------------------------------------------------------
    
PLOT_DF <- left_join(MAGMA_DF, LDSR_DF,
                         by = 'Category') %>% reshape2::melt() %>%
  mutate(cell_type = case_when(Category %in% c("CCK_RELN", "CCK_SORCS1", "CCK_SYT6", "ID2_CSMD1", "LAMP5_CCK", "LAMP5_NDNF", "LAMP5_NOS1", "PV_SCUBE3", "PV_SST", "PV_SULF1", "PV_WFDC2", "SST_ADGRG6", "SST_B3GAT2", "SST_BRINP3", "SST_CALB1", "SST_NPY", "SST_STK32A", "SST_TH", "VIP_ABI3BP", "VIP_ADAMTSL1", "VIP_CHRM2", "VIP_CRH", "VIP_DPP6", "VIP_HS3ST3A1", "VIP_KIRREL3", "VIP_PCDH20") ~ "Mature Inhibitory Neurons",
                                   Category %in% c("CGE_dev", "ID2_dev", "MGE_dev-1", "MGE_dev-2", "PV_dev", "PV_SCUBE3_dev", "PV_SULF1_dev", "SST_ADGRG6_dev", "SST_CALB1_dev", "VIP_dev") ~ "Developing Inhibitory Neurons",
                                   Category %in% c("BKGR_NRGN", "L2_CUX2_LAMP5", "L3_CUX2_PRSS12", "L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_MME", "L5-6_THEMIS_CNR1", "L5-6_THEMIS_NTNG2", "L5-6_TLE4_HTR2C", "L5-6_TLE4_SCUBE1", "L5-6_TLE4_SORCS1") ~ "Mature Excitatory Neurons",
                                   Category %in% c("L2-3_CUX2_dev-1", "L2-3_CUX2_dev-2", "L2-3_CUX2_dev-3", "L2-3_CUX2_dev-4", "L2-3_CUX2_dev-5", "L2-3_CUX2_dev-6", "L2-3_CUX2_dev-fetal", "L2_CUX2_LAMP5_dev", "L4_RORB_dev-1", "L4_RORB_dev-2", "L4_RORB_dev-fetal", "L5-6_THEMIS_dev-1", "L5-6_THEMIS_dev-2", "L5-6_TLE4_dev", "PN_dev") ~ "Developing Excitatory Neurons",
                                   Category %in% c("Astro_dev-1", "Astro_dev-2", "Astro_dev-3", "Astro_dev-4", "Astro_dev-5", "Astro_GFAP", "Astro_SLC1A2", "Astro_SLC1A2_dev", "Micro", "Micro_out", "Oligo-1", "Oligo-2", "Oligo-3", "Oligo-4", "Oligo-5", "Oligo_mat", "OPC", "OPC_dev", "OPC_MBP", "Vas_CLDN5", "Vas_PDGFRB", "Vas_TBX18") ~ "Non-Neuronal")) %>%
  mutate(cell_type=factor(cell_type, levels = c("Mature Excitatory Neurons", "Developing Excitatory Neurons", "Mature Inhibitory Neurons", "Developing Inhibitory Neurons", "Non-Neuronal")))
    
PLOT_DF$Category <- gsub("_", "-", PLOT_DF$Category)
    
MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))), width = 0.8,
                                                  fill = variable, group = rev(variable))) +
      geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(GWAS) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 16),
            axis.title.y = element_text(colour = "#000000", size = 16),
            axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 13),
            legend.title = element_blank(),
            title = element_text(colour = "#000000", size = 16),
            legend.text = element_text(colour = "#000000", size = 15),
            strip.text = element_text(size=14, face = 'bold')) + #,
      #legend.position = "top") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 8) +
      facet_col(facets = vars(cell_type), scales = "free", space = "free", strip.position = "top")
    
    assign(paste0(GWAS, '_magma_ldsr_', PROJECT, '_lvl', level, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv)
    
  }
}

## Save plots -------------------------------------------------------------------------

legend <- get_legend(HEIGHT_magma_ldsr_herring_lvl2_plot)

#Title plots
NEUROTICISM_magma_ldsr_herring_lvl2_plot <- NEUROTICISM_magma_ldsr_herring_lvl2_plot + ggtitle("Neuroticism") 
ASD_magma_ldsr_herring_lvl2_plot <- ASD_magma_ldsr_herring_lvl2_plot + ggtitle("Autism") 
HEIGHT_magma_ldsr_herring_lvl2_plot <- HEIGHT_magma_ldsr_herring_lvl2_plot + ggtitle("Height") 

#Save plots

jpeg(file = paste0(FIG_DIR, 'HEIGHT_magma_ldsr_herring_plot.jpeg'), units = "in", width = 8.75, height = 22, res = 300)
plot(HEIGHT_magma_ldsr_herring_lvl2_plot)
dev.off() 

jpeg(file = paste0(FIG_DIR, 'BPD_magma_ldsr_herring_plot.jpeg'), units = "in", width = 8.75, height = 22, res = 300)
plot(BPD_magma_ldsr_herring_lvl2_plot)
dev.off() 

jpeg(file = paste0(FIG_DIR, 'ASD_magma_ldsr_herring_plot.jpeg'), units = "in", width = 8.75, height = 22, res = 300)
plot(ASD_magma_ldsr_herring_lvl2_plot)
dev.off() 

jpeg(file = paste0(FIG_DIR, 'ADHD_magma_ldsr_herring_plot.jpeg'), units = "in", width = 8.75, height = 22, res = 300)
plot(ADHD_magma_ldsr_herring_lvl2_plot)
dev.off() 

jpeg(file = paste0(FIG_DIR, 'MDD_magma_ldsr_herring_plot.jpeg'), units = "in", width = 8.75, height = 22, res = 300)
plot(MDD_magma_ldsr_herring_lvl2_plot)
dev.off() 

jpeg(file = paste0(FIG_DIR, 'NEUROTICISM_magma_ldsr_herring_plot.jpeg'), units = "in", width = 8.75, height = 22, res = 300)
plot(NEUROTICISM_magma_ldsr_herring_lvl2_plot)
dev.off()
